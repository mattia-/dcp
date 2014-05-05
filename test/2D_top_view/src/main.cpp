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
    parameters.add ("target_u_file_name", "u_target");
    parameters.add ("target_p_file_name", "p_target");
    
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
    dolfin::Function targetU (completeV[0] -> collapse ());
    dolfin::Function targetP (completeV[1] -> collapse ());
    
    dolfin::HDF5File targetUFile (static_cast<std::string> (parameters ["target_u_file_name"]) + ".hdf5", "r");
    targetUFile.read (targetU, "u");
    dolfin::HDF5File targetPFile (static_cast<std::string> (parameters ["target_p_file_name"]) + ".hdf5", "r");
    targetPFile.read (targetP, "p");
    
    
    // create mesh and finite element space on partial geometry
    dolfin::Mesh mesh (parameters ["partial_mesh_file_name"]);
    primal::FunctionSpace V (mesh);
    
    // interpolate target functions from completeMesh to actual mesh
    dolfin::Function interpolatedTargetU (V[0] -> collapse ());
    dolfin::Function interpolatedTargetP (V[1] -> collapse ());
    
    interpolatedTargetU.interpolate (targetU);
    interpolatedTargetP.interpolate (targetP);

    
    // ============================================================================ //
    // =========================== COMPOSITE PROBLEM ============================== //
    // ============================================================================ //
    // define problems
    DCP::NonlinearDifferentialProblem <primal::ResidualForm, primal::JacobianForm> 
        primalProblem (dolfin::reference_to_no_delete_pointer (mesh), 
                        dolfin::reference_to_no_delete_pointer (V),
                        "trial");

    DCP::LinearDifferentialProblem <adjoint::BilinearForm, adjoint::LinearForm> 
        adjointProblem (dolfin::reference_to_no_delete_pointer (mesh),
                         dolfin::reference_to_no_delete_pointer (V));

    
    // create composite differential problem
    DCP::CompositeDifferentialProblem problems;
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
    problems["adjoint"].setCoefficient ("linear_form", dolfin::reference_to_no_delete_pointer (targetU), "U");
    problems["adjoint"].setCoefficient ("linear_form", dolfin::reference_to_no_delete_pointer (targetP), "P");

    problems["adjoint"].addDirichletBC (dolfin::DirichletBC (*V[0], adjoint_dirichletBC, adjoint_dirichletBoundary));

    problems["adjoint"].setIntegrationSubdomains ("bilinear_form", 
                                                  dolfin::reference_to_no_delete_pointer (meshFacets),
                                                  DCP::SubdomainType::BOUNDARY_FACETS);
    problems["adjoint"].setIntegrationSubdomains ("linear_form", 
                                                  dolfin::reference_to_no_delete_pointer (meshCells),
                                                  DCP::SubdomainType::INTERNAL_CELLS);

    problems.addLink ("adjoint", "u", "bilinear_form", "primal", 0);
    problems.addLink ("adjoint", "u", "linear_form", "primal", 0);
    problems.addLink ("adjoint", "p", "linear_form", "primal", 1);
    

    // ============================================================================== //
    // =========================== OBJECTIVE FUNCTIONAL ============================= //
    // ============================================================================== //
    // define functional
    DCP::ObjectiveFunctional <objective_functional::Form_J, objective_functional::Gradient>
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

    // define subdomains for objective functional
    objective_functional::ControlDomain objective_functional_controlDomain;
    objective_functional_controlDomain.mark (meshFacets, 2); 
    
    // functional settings
    objectiveFunctional.setCoefficient ("functional", dolfin::reference_to_no_delete_pointer (sigma_1), "sigma_1");
    objectiveFunctional.setCoefficient ("functional", dolfin::reference_to_no_delete_pointer (sigma_2), "sigma_2");
    objectiveFunctional.setCoefficient ("functional", dolfin::reference_to_no_delete_pointer (g), "g");
    objectiveFunctional.setCoefficient ("functional", dolfin::reference_to_no_delete_pointer (targetU), "U");
    objectiveFunctional.setCoefficient ("functional", dolfin::reference_to_no_delete_pointer (targetP), "P");
    objectiveFunctional.setCoefficient ("functional", 
                                        dolfin::reference_to_no_delete_pointer (problems["primal"].solution()[0]), 
                                        "u");
    objectiveFunctional.setCoefficient ("functional", 
                                        dolfin::reference_to_no_delete_pointer (problems["primal"].solution()[1]), 
                                        "p");
    
    objectiveFunctional.setIntegrationSubdomains (dolfin::reference_to_no_delete_pointer (meshFacets), 
                                                  DCP::SubdomainType::BOUNDARY_FACETS);
    objectiveFunctional.setIntegrationSubdomains (dolfin::reference_to_no_delete_pointer (meshCells),
                                                  DCP::SubdomainType::INTERNAL_CELLS);
    
    objectiveFunctional.setCoefficient ("gradient",
                                        dolfin::reference_to_no_delete_pointer (problems["adjoint"].solution()[1]),
                                        "theta");
    objectiveFunctional.setCoefficient ("gradient",
                                        dolfin::reference_to_no_delete_pointer (g),
                                        "g");

    
    
    // ============================================================================== //
    // ===================== OBJECTIVE FUNCTIONAL COMPONENTS ======================== //
    // ============================================================================== //
    // ----------------------
    // component: only target
    // ----------------------
    objective_functional::Form_J_control_component functionalControlComponent (mesh);
        
    // functional settings
    functionalControlComponent.set_coefficient ("sigma_1", dolfin::reference_to_no_delete_pointer (sigma_1));
    functionalControlComponent.set_coefficient ("sigma_2", dolfin::reference_to_no_delete_pointer (sigma_2));
    functionalControlComponent.set_coefficient ("g", dolfin::reference_to_no_delete_pointer (g));
    
    functionalControlComponent.set_exterior_facet_domains (dolfin::reference_to_no_delete_pointer (meshFacets));
    
    
    // -----------------------
    // component: only control
    // -----------------------
    objective_functional::Form_J_target_component functionalTargetComponent (mesh);
        
    // functional settings
    functionalTargetComponent.set_coefficient ("U", dolfin::reference_to_no_delete_pointer (targetU));
    functionalTargetComponent.set_coefficient ("P", dolfin::reference_to_no_delete_pointer (targetP));
    functionalTargetComponent.set_coefficient ("u", dolfin::reference_to_no_delete_pointer (problems["primal"].solution()[0]));
    functionalTargetComponent.set_coefficient ("p", dolfin::reference_to_no_delete_pointer (problems["primal"].solution()[1]));
    
    functionalTargetComponent.set_cell_domains (dolfin::reference_to_no_delete_pointer (meshCells));
    
    

    // ============================================================================== //
    // =============================== LIFT AND DRAG ================================ //
    // ============================================================================== //
    primal_noSlipBoundary.mark (meshFacets, 3);
    
    // ------
    // lift
    // ------
    lift_drag::Form_Lift liftComputer (mesh);
        
    // functional settings
    liftComputer.set_coefficient ("p", dolfin::reference_to_no_delete_pointer (problems["primal"].solution()[1]));
    liftComputer.set_exterior_facet_domains (dolfin::reference_to_no_delete_pointer (meshFacets));
    
    
    // ------
    // drag
    // ------
    lift_drag::Form_Drag dragComputer (mesh);
        
    // functional settings
    dragComputer.set_coefficient ("p", dolfin::reference_to_no_delete_pointer (problems["primal"].solution()[1]));
    dragComputer.set_exterior_facet_domains (dolfin::reference_to_no_delete_pointer (meshFacets));
    
    

    // ============================================================================== //
    // =============================== OPTIMIZATION  ================================ //
    // ============================================================================== //
    ControlDirichletBC controlDirichletBC (g);
    problems["primal"].addDirichletBC (dolfin::DirichletBC (*(*V[0])[0], controlDirichletBC, primal_inflowBoundary), "x_inflow_BC");
    
    // define control value updater
    ValueUpdater updater;
    
    // define search direction computer
    objective_functional::SearchDirectionComputer searchDirectionComputer (controlMesh);
    
    // define optimizer
    DCP::BacktrackingOptimizer backtrackingOptimizer;
    backtrackingOptimizer.parameters ["relative_increment_tolerance"] = 1e-5;
//    backtrackingOptimizer.parameters ["max_minimization_iterations"] = 200;
    backtrackingOptimizer.parameters ["output_file_name"] = "results.txt";
    
    backtrackingOptimizer.apply (problems, objectiveFunctional, g, updater, searchDirectionComputer);
            

    // solve problem with computed control value
    updater (problems, g);
    
    problems.solve ();
    
    dolfin::plot (problems.solution ("primal")[0], "solution of the problem with computed control. Velocity");
    dolfin::plot (problems.solution ("primal")[0][0], "solution of the problem with computed control. x velocity");
    dolfin::plot (problems.solution ("primal")[0][1], "solution of the problem with computed control. y velocity");
    dolfin::plot (problems.solution ("primal")[1], "solution of the problem with computed control. Pressure");
    dolfin::plot (problems.solution ("adjoint")[1], "solution of the adjoint problem with computed control. Pressure");
    dolfin::plot (g, "Control");
    
    dolfin::Function differenceU (interpolatedTargetU);
    dolfin::Function differenceU1 (interpolatedTargetU);
    differenceU1.interpolate (problems.solution ("primal")[0]);
    differenceU = differenceU1 - interpolatedTargetU; 
    dolfin::Function differenceP (interpolatedTargetP);
    dolfin::Function differenceP1 (interpolatedTargetP);
    differenceP1.interpolate (problems.solution ("primal")[1]);
    differenceP = differenceP1 - interpolatedTargetP;
    
    
    dolfin::plot (meshCells, "Control region");
    dolfin::plot (differenceU, "Difference between target and controlled velocity");
    dolfin::plot (differenceP, "Difference between target and controlled pressure");
    
    dolfin::interactive ();
    
    
    // compute functional components
    std::ofstream componentsOutputStream ("functional_components.txt");
    componentsOutputStream << "Control_component = " << dolfin::assemble (functionalControlComponent) << std::endl;
    componentsOutputStream << "Target_component = " << dolfin::assemble (functionalTargetComponent) << std::endl;
    componentsOutputStream.close ();
    

    // compute lift and drag
    std::ofstream liftDragOutputStream ("reconstructed_lift_drag.txt");
    liftDragOutputStream << "Lift = " << dolfin::assemble (liftComputer) << std::endl;
    liftDragOutputStream << "Drag = " << dolfin::assemble (dragComputer) << std::endl;
    liftDragOutputStream.close ();
    return 0;
}
