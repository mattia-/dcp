/* 
 *  Copyright (C) 2014, Mattia Tamellini, mattia.tamellini@gmail.com
 * 
 *  This file is part of the DCP library
 *   
 *   The DCP library is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   The DCP library is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with the DCP library.  If not, see <http://www.gnu.org/licenses/>. 
 */ 

#include <iostream>
#include <string>
#include <dolfin.h>
#include <DifferentialProblem/DifferentialProblem.hpp>
#include <ObjectiveFunctional/Functional.hpp>
#include <Optimizer/Optimizer.hpp>
#include "primal.h"
#include "adjoint.h"
#include "objective_functional.h"
#include "control_variable_function_space.h"
#include "lift_drag.h"
#include "main_settings.h"

int main (int argc, char* argv[])
{
    // ============================================================================ //
    // =========================== PROBLEM PARAMETERS ============================= //
    // ============================================================================ //
    // log level
//    dolfin::set_log_level (dolfin::DBG);
//    dolfin::set_log_level (dolfin::PROGRESS);
//    dolfin::set_log_level (dolfin::WARNING);

    // define parameters and their default values
    dolfin::Parameters parameters ("main_parameters");
    parameters.add ("complete_mesh_file_name", "../src/complete_mesh/complete_mesh.xml");
    parameters.add ("partial_mesh_file_name", "../src/partial_mesh/partial_mesh.xml");
    parameters.add ("target_solution_file_name", "target_solution");
    
    // the next parameter is set to allow projection of solution from complete mesh to partial mesh.
    // In any case, the values that will be extrapolated are not important because the functional will
    // integrate the difference of the two only on a subdomain, which does not include the points whre
    // extrapolation takes place
    dolfin::parameters ["allow_extrapolation"] = true;
    
    // read parameters from command line and overwrite default values
    parameters.parse (argc, argv);
    
    // ============================================================================ //
    // =========================== TARGET SOLUTION ================================ //
    // ============================================================================ //
    // create mesh and finite element space to read target solution
    dolfin::Mesh completeMesh (parameters ["complete_mesh_file_name"]);
    primal::FunctionSpace completeV (completeMesh);
    
    // read target solutions from file
    dolfin::Function targetSolution (completeV);
    
    dolfin::HDF5File targetSolutionFile (MPI_COMM_WORLD, 
                                         static_cast<std::string> (parameters ["target_solution_file_name"]) + ".hdf5", 
                                         "r");
    targetSolutionFile.read (targetSolution, "solution");
    
    
    // create mesh and finite element space on partial geometry
    dolfin::Mesh mesh (parameters ["partial_mesh_file_name"]);
    primal::FunctionSpace V (mesh);
    
    // interpolate target functions from completeMesh to actual mesh
    dolfin::Function interpolatedTargetSolution (V);
    interpolatedTargetSolution.interpolate (targetSolution);

    
    // ============================================================================ //
    // =========================== COMPOSITE PROBLEM ============================== //
    // ============================================================================ //
    // define problems
    dcp::NonlinearDifferentialProblem <primal::ResidualForm, primal::JacobianForm> 
        primalProblem (dolfin::reference_to_no_delete_pointer (mesh), 
                        dolfin::reference_to_no_delete_pointer (V),
                        "trial");

    dcp::LinearDifferentialProblem <adjoint::BilinearForm, adjoint::LinearForm> 
        adjointProblem (dolfin::reference_to_no_delete_pointer (mesh),
                         dolfin::reference_to_no_delete_pointer (V));

    
    // create composite differential problem
    dcp::CompositeDifferentialProblem problems;
    problems.addProblem ("primal", primalProblem);
    problems.addProblem ("adjoint", adjointProblem);
    
    // define constants
    dolfin::Constant nu (1e-1);
    dolfin::Constant primal_yInflowDirichletBC (0.0);
    dolfin::Constant primal_symmetryDirichletBC (0.0);
    dolfin::Constant primal_noSlipCondition (0.0, 0.0);
    dolfin::Constant adjoint_dirichletBC (0.0, 0.0);

    // define boundary conditions subdomains
    primal::InflowBoundary primal_inflowBoundary;
    primal::GammaSD primal_gammaSD;
    primal::NoSlipBoundary primal_noSlipBoundary;
    adjoint::DirichletBoundary adjoint_dirichletBoundary;
    adjoint::RobinBoundary adjoint_robinBoundary;
    adjoint::ExternalLoadDomain adjoint_externalLoadDomain;
    
    // define intergration subdomains
    dolfin::FacetFunction<std::size_t> meshFacets (mesh);
    meshFacets.set_all (0);
    adjoint_robinBoundary.mark (meshFacets, 1);
    
    dolfin::CellFunction<std::size_t> meshCells (mesh);
    meshCells.set_all (0);
    adjoint_externalLoadDomain.mark (meshCells, 1);

    
    // primal problem settings
    problems["primal"].setCoefficient ("residual_form", dolfin::reference_to_no_delete_pointer (nu), "nu");
    problems["primal"].setCoefficient ("jacobian_form", dolfin::reference_to_no_delete_pointer (nu), "nu");

    problems["primal"].addDirichletBC (dolfin::DirichletBC (*(*V[0])[1], primal_yInflowDirichletBC, primal_inflowBoundary));
    problems["primal"].addDirichletBC (dolfin::DirichletBC (*V[0], primal_noSlipCondition, primal_noSlipBoundary));
    problems["primal"].addDirichletBC (dolfin::DirichletBC (*(*V[0])[1], primal_symmetryDirichletBC, primal_gammaSD));


    // adjoint problem settings
    problems["adjoint"].setCoefficient ("bilinear_form", dolfin::reference_to_no_delete_pointer (nu), "nu");
    problems["adjoint"].setCoefficient ("linear_form", dolfin::reference_to_no_delete_pointer (targetSolution [0]), "U");
    problems["adjoint"].setCoefficient ("linear_form", dolfin::reference_to_no_delete_pointer (targetSolution [1]), "P");

    problems["adjoint"].addDirichletBC (dolfin::DirichletBC (*V[0], adjoint_dirichletBC, adjoint_dirichletBoundary));

    problems["adjoint"].setIntegrationSubdomains ("bilinear_form", 
                                                  dolfin::reference_to_no_delete_pointer (meshFacets),
                                                  dcp::SubdomainType::BOUNDARY_FACETS);
    problems["adjoint"].setIntegrationSubdomains ("linear_form", 
                                                  dolfin::reference_to_no_delete_pointer (meshCells),
                                                  dcp::SubdomainType::INTERNAL_CELLS);

    problems.addLink ("adjoint", "u", "bilinear_form", "primal", 0);
    problems.addLink ("adjoint", "u", "linear_form", "primal", 0);
    problems.addLink ("adjoint", "p", "linear_form", "primal", 1);
    

    // ============================================================================== //
    // =========================== OBJECTIVE FUNCTIONAL ============================= //
    // ============================================================================== //
    // define functional
    dcp::ObjectiveFunctional <objective_functional::Form_J, objective_functional::Gradient>
        objectiveFunctional (mesh);
    
    // define functional coefficients
    dolfin::Constant sigma_1 (objective_functional::sigma_1);
    dolfin::Constant sigma_2 (objective_functional::sigma_2);
    
    // define control variable
    // mesh:
    dolfin::IntervalMesh controlMesh (100, 0, 7);
    // function space:
    control_variable_function_space::FunctionSpace controlFunctionSpace (controlMesh);
    // control variable itself:
    dolfin::Function g (controlFunctionSpace);
    g = dolfin::Constant (1.0);
    ControlDirichletBC controlDirichletBC (g);

    // define subdomains for objective functional
    objective_functional::ControlDomain objective_functional_controlDomain;
    objective_functional_controlDomain.mark (meshFacets, 2); 
    
    // functional settings
    objectiveFunctional.setCoefficient ("functional", dolfin::reference_to_no_delete_pointer (sigma_1), "sigma_1");
    objectiveFunctional.setCoefficient ("functional", dolfin::reference_to_no_delete_pointer (sigma_2), "sigma_2");
    objectiveFunctional.setCoefficient ("functional", dolfin::reference_to_no_delete_pointer (controlDirichletBC), "g");
    objectiveFunctional.setCoefficient ("functional", dolfin::reference_to_no_delete_pointer (targetSolution [0]), "U");
    objectiveFunctional.setCoefficient ("functional", dolfin::reference_to_no_delete_pointer (targetSolution [1]), "P");
    objectiveFunctional.setCoefficient ("functional", 
                                        dolfin::reference_to_no_delete_pointer (problems.solution("primal")[0]), 
                                        "u");
    objectiveFunctional.setCoefficient ("functional", 
                                        dolfin::reference_to_no_delete_pointer (problems.solution("primal")[1]), 
                                        "p");
    
    objectiveFunctional.setIntegrationSubdomains (dolfin::reference_to_no_delete_pointer (meshFacets), 
                                                  dcp::SubdomainType::BOUNDARY_FACETS);
    objectiveFunctional.setIntegrationSubdomains (dolfin::reference_to_no_delete_pointer (meshCells),
                                                  dcp::SubdomainType::INTERNAL_CELLS);
    
    objectiveFunctional.setCoefficient ("gradient",
                                        dolfin::reference_to_no_delete_pointer (problems.solution("adjoint")[1]),
                                        "theta");
    objectiveFunctional.setCoefficient ("gradient",
                                        dolfin::reference_to_no_delete_pointer (controlDirichletBC),
                                        "g");

    
    
    // ============================================================================== //
    // ===================== OBJECTIVE FUNCTIONAL COMPONENTS ======================== //
    // ============================================================================== //
    // ----------------------
    // component: velocity
    // ----------------------
    objective_functional::Form_J_u functionalVelocityComponent (mesh);
        
    // functional settings
    functionalVelocityComponent.set_coefficient ("U", dolfin::reference_to_no_delete_pointer (targetSolution [0]));
    functionalVelocityComponent.set_coefficient ("u", 
                                                 dolfin::reference_to_no_delete_pointer (problems.solution("primal")[0]));
    
    functionalVelocityComponent.set_cell_domains (dolfin::reference_to_no_delete_pointer (meshCells));
    
    // ----------------------
    // component: pressure
    // ----------------------
    objective_functional::Form_J_p functionalPressureComponent (mesh);
        
    // functional settings
    functionalPressureComponent.set_coefficient ("P", dolfin::reference_to_no_delete_pointer (targetSolution [1]));
    functionalPressureComponent.set_coefficient ("p", dolfin::reference_to_no_delete_pointer (problems.solution("primal")[1]));
    
    functionalPressureComponent.set_cell_domains (dolfin::reference_to_no_delete_pointer (meshCells));
    
    // ----------------------
    // component: control
    // ----------------------
    objective_functional::Form_J_g functionalControlComponent (mesh);
        
    // functional settings
    functionalControlComponent.set_coefficient ("sigma_1", dolfin::reference_to_no_delete_pointer (sigma_1));
    functionalControlComponent.set_coefficient ("g", dolfin::reference_to_no_delete_pointer (controlDirichletBC));
    
    functionalControlComponent.set_exterior_facet_domains (dolfin::reference_to_no_delete_pointer (meshFacets));
    
    // -----------------------------
    // component: control derivative
    // -----------------------------
    objective_functional::Form_J_gDer functionalControlDerivativeComponent (mesh);
        
    // functional settings
    functionalControlDerivativeComponent.set_coefficient ("sigma_2", dolfin::reference_to_no_delete_pointer (sigma_2));
    functionalControlDerivativeComponent.set_coefficient ("g", dolfin::reference_to_no_delete_pointer (controlDirichletBC));
    
    functionalControlDerivativeComponent.set_exterior_facet_domains (dolfin::reference_to_no_delete_pointer (meshFacets));
    
    
    // compute and print to file functional components
    std::ofstream componentsOutputStream ("functional_components.txt");
    componentsOutputStream << "INITIAL VALUES:" << std::endl;
    componentsOutputStream << "Velocity_component = " << dolfin::assemble (functionalVelocityComponent) << std::endl;
    componentsOutputStream << "Pressure_component = " << dolfin::assemble (functionalPressureComponent) << std::endl;
    componentsOutputStream << "Control_component = " << dolfin::assemble (functionalControlComponent) << std::endl;
    componentsOutputStream << "Control_derivative_component = " 
                           << dolfin::assemble (functionalControlDerivativeComponent) << std::endl;
    componentsOutputStream << std::endl;

    
    // ============================================================================== //
    // =============================== LIFT AND DRAG ================================ //
    // ============================================================================== //
    primal_noSlipBoundary.mark (meshFacets, 3);
    
    // ------
    // lift
    // ------
    lift_drag::Form_Lift liftComputer (mesh);
        
    // functional settings
    liftComputer.set_coefficient ("p", dolfin::reference_to_no_delete_pointer (problems.solution("primal")[1]));
    liftComputer.set_exterior_facet_domains (dolfin::reference_to_no_delete_pointer (meshFacets));
    
    
    // ------
    // drag
    // ------
    lift_drag::Form_Drag dragComputer (mesh);
        
    // functional settings
    dragComputer.set_coefficient ("p", dolfin::reference_to_no_delete_pointer (problems.solution("primal")[1]));
    dragComputer.set_exterior_facet_domains (dolfin::reference_to_no_delete_pointer (meshFacets));
    

    // ============================================================================== //
    // =============================== OPTIMIZATION  ================================ //
    // ============================================================================== //
    problems["primal"].addDirichletBC (dolfin::DirichletBC (*(*V[0])[0], controlDirichletBC, primal_inflowBoundary), "x_inflow_BC");
    
    // define control value updater
    ValueUpdater updater;
    
    // define search direction computer
    objective_functional::SearchDirectionComputer searchDirectionComputer (controlMesh);
    
    // define optimizer
    dcp::BacktrackingOptimizer backtrackingOptimizer;
    backtrackingOptimizer.parameters ["output_file_name"] = "results.txt";
    backtrackingOptimizer.parameters ["relative_increment_tolerance"] = 1e-3;
    
    backtrackingOptimizer.apply (problems, objectiveFunctional, g, updater, searchDirectionComputer);
            

    // compute difference between target and reconstructed solution
    dolfin::Function difference (V);
    difference = problems.solution ("primal") - interpolatedTargetSolution; 
    
    
    // ============================================================================== //
    // =============================== POST PROCESSING ============================== //
    // ============================================================================== //
    // plots
    dolfin::VTKPlotter meshPlotter (dolfin::reference_to_no_delete_pointer (mesh));
    meshPlotter.parameters["title"] = "Mesh";
    meshPlotter.plot ();
    
    dolfin::VTKPlotter primalVelocityPlotter (dolfin::reference_to_no_delete_pointer (problems.solution ("primal")[0]));
    primalVelocityPlotter.parameters["title"] = "Solution of the problem with computed control. Velocity";
    primalVelocityPlotter.parameters["input_keys"] = "m";
    primalVelocityPlotter.plot ();
    
    dolfin::VTKPlotter primalPressurePlotter (dolfin::reference_to_no_delete_pointer (problems.solution ("primal")[1]));
    primalPressurePlotter.parameters["mode"]  = "color";
    primalPressurePlotter.parameters["title"] = "Solution of the problem with computed control. Pressure";
    primalPressurePlotter.plot ();
    
    dolfin::VTKPlotter adjointPressurePlotter (dolfin::reference_to_no_delete_pointer (problems.solution ("adjoint")[1]));
    adjointPressurePlotter.parameters["mode"]  = "color";
    adjointPressurePlotter.parameters["title"] = "Solution of the adjoint problem with computed control. Pressure";
    adjointPressurePlotter.plot ();
    
    dolfin::VTKPlotter controlPlotter (dolfin::reference_to_no_delete_pointer (g));
    controlPlotter.parameters["title"] = "Control";
    controlPlotter.plot ();
    
    dolfin::VTKPlotter controlRegionPlotter (dolfin::reference_to_no_delete_pointer (meshCells));
    controlRegionPlotter.parameters["mode"]  = "color";
    controlRegionPlotter.parameters["title"] = "Control region";
    controlRegionPlotter.plot ();
    
    dolfin::VTKPlotter velocityDifferencePlotter (dolfin::reference_to_no_delete_pointer (difference [0]));
    velocityDifferencePlotter.parameters["title"] = "Difference between target and reconstructed velocity";
    velocityDifferencePlotter.parameters["mode"]  = "color";
    velocityDifferencePlotter.plot ();
    
    dolfin::VTKPlotter pressureDifferencePlotter (dolfin::reference_to_no_delete_pointer (difference [1]));
    pressureDifferencePlotter.parameters["mode"]  = "color";
    pressureDifferencePlotter.parameters["title"] = "Difference between target and reconstructed pressure";
    pressureDifferencePlotter.plot ();
    
    dolfin::interactive ();
    
    
    // print results to file
    dolfin::File velocityOutputFile ("reconstructed_velocity.pvd");
    velocityOutputFile << problems.solution ("primal")[0];
    
    dolfin::File pressureOutputFile ("reconstructed_pressure.pvd");
    pressureOutputFile << problems.solution ("primal")[1];
    
    dolfin::File velocityDifferenceFile ("velocity_difference.pvd");
    velocityDifferenceFile << difference [0];
    
    dolfin::File pressureDifferenceFile ("pressure_difference.pvd");
    pressureDifferenceFile << difference [1];

    
    // print control variable to file
    std::ofstream gFile ("control_values.txt");
    for (dolfin::CellIterator c(controlMesh); !c.end(); ++c)
    {
        for (dolfin::VertexIterator v0(*c); !v0.end(); ++v0)
        {
            dolfin::Array<double> x (1);
            dolfin::Array<double> gValues (1);
            x[0] = v0 -> x(0);
            g.eval (gValues, x);
            gFile << x [0] << " " << gValues [0] << std::endl;
        }
    }
    
    
    // compute and print to file functional components
    componentsOutputStream << "FINAL VALUES:" << std::endl;
    componentsOutputStream << "Velocity_component = " << dolfin::assemble (functionalVelocityComponent) << std::endl;
    componentsOutputStream << "Pressure_component = " << dolfin::assemble (functionalPressureComponent) << std::endl;
    componentsOutputStream << "Control_component = " << dolfin::assemble (functionalControlComponent) << std::endl;
    componentsOutputStream << "Control_derivative_component = " 
                           << dolfin::assemble (functionalControlDerivativeComponent) << std::endl;
    componentsOutputStream.close ();
    

    // compute and print to file lift and drag
    std::ofstream liftDragOutputStream ("reconstructed_lift_drag.txt");
    liftDragOutputStream << "Lift = " << dolfin::assemble (liftComputer) << std::endl;
    liftDragOutputStream << "Drag = " << dolfin::assemble (dragComputer) << std::endl;
    liftDragOutputStream.close ();
    return 0;
}
