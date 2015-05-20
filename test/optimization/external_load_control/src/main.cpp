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
#include <mshr.h>
#include <dcp/differential_problems/differential_problems.h>
#include <dcp/objective_functional/objective_functional.h>
#include <dcp/expressions/expressions.h>
#include <dcp/optimizers/optimizers.h>
#include "primal.h"
#include "adjoint.h"
#include "objective_functional.h"

namespace objective_functional
{
    class Gradient : public dcp::VariableExpression
    {
        void eval (dolfin::Array<double>& values, const dolfin::Array<double>& x) const
        {
            evaluateVariable ("p", values, x);
            values[0] = -values[0];
            values[1] = -values[1]; // TODO check gradient expression
        }
    };
}

class ProblemUpdater
{
    public:
        void operator() (dcp::EquationSystem& problem, 
                         const dolfin::GenericFunction& controlValue) const
        {
            problem["primal"].setCoefficient ("linear_form", 
                                              dolfin::reference_to_no_delete_pointer (controlValue),
                                              "g");
        }
};

class DirichletBoundary
{
    public:
        bool operator() (const dolfin::Array<double>& x, bool on_boundary)
        {
            return on_boundary;
        }
};

class ObservationDomain
{
    public:
        bool operator() (const dolfin::Array<double>& x, bool on_boundary)
        {
            return x[0] <= 0.75 && x[0] >= 0.25 && x[1] <= 0.75 && x[1] >= 0.25;
        }
};

    
int main (int argc, char* argv[])
{
    // log level
    // dolfin::set_log_level (dolfin::DBG);
    
    // mesh
    dolfin::UnitSquareMesh mesh (20, 20);
    
    // function spaces
    primal::FunctionSpace V (mesh);
    adjoint::FunctionSpace W (mesh);

    
    // =======================
    // PRIMAL+ADJOINT PROBLEMS
    // =======================
    // define problems
    dcp::LinearProblem <primal::BilinearForm, primal::LinearForm> 
        primalProblem (dolfin::reference_to_no_delete_pointer (V));

    dcp::LinearProblem <adjoint::BilinearForm, adjoint::LinearForm> 
        adjointProblem (dolfin::reference_to_no_delete_pointer (W));

    // create composite differential problem
    dcp::EquationSystem problems;
    problems.addProblem ("primal", primalProblem);
    problems.addProblem ("adjoint", adjointProblem);
    
    // define constants
    dolfin::Constant u_0 (1.0);
    dolfin::Constant dirichletBC (0.0);

    // define boundary conditions subdomains
    dcp::Subdomain dirichletBoundary ((DirichletBoundary ()));
    
    // define observation subdomain
    dcp::Subdomain observationDomain ((ObservationDomain ()));
    
    // define intergration subdomains
    dolfin::CellFunction<std::size_t> meshCells (mesh);
    meshCells.set_all (0);
    observationDomain.mark (meshCells, 1);
    
    // primal problem settings
    problems["primal"].addDirichletBC (dirichletBC, dirichletBoundary);

    // adjoint problem settings
    problems["adjoint"].addDirichletBC (dirichletBC, dirichletBoundary);
    problems["adjoint"].setCoefficient ("linear_form", dolfin::reference_to_no_delete_pointer (u_0), "u_0");
    problems["adjoint"].setIntegrationSubdomain ("linear_form", 
                                                  dolfin::reference_to_no_delete_pointer (meshCells),
                                                  dcp::SubdomainType::INTERNAL_CELLS);

    problems.addLink ("adjoint", "u", "linear_form", "primal");
    

    // ====================
    // OBJECTIVE FUNCTIONAL 
    // ====================
    // define functional
    dcp::ObjectiveFunctional <objective_functional::Form_J, objective_functional::Gradient>
        objectiveFunctional (mesh);
    
    // control variable 
    dolfin::Function g (V);
    g = dolfin::Constant (0.0);
    
    // functional settings
    objectiveFunctional.setCoefficient ("functional", dolfin::reference_to_no_delete_pointer (u_0), "u_0");
    objectiveFunctional.setCoefficient ("functional", 
                                        dolfin::reference_to_no_delete_pointer (problems.solution ("primal")), 
                                        "u");
    
    objectiveFunctional.setIntegrationSubdomain (dolfin::reference_to_no_delete_pointer (meshCells),
                                                  dcp::SubdomainType::INTERNAL_CELLS);
    
    objectiveFunctional.setCoefficient ("gradient",
                                        dolfin::reference_to_no_delete_pointer (problems.solution ("adjoint")),
                                        "p");
    
    
    // ============
    // OPTIMIZATION
    // ============
    problems["primal"].setCoefficient ("linear_form", dolfin::reference_to_no_delete_pointer (g), "g");
    
    // define control value updater
    ProblemUpdater updater;
    
    // define optimizer
    dcp::BacktrackingOptimizer backtrackingOptimizer;
    backtrackingOptimizer.parameters ["output_file_name"] = "results.txt";
    backtrackingOptimizer.parameters ["relative_increment_tolerance"] = 1e-3;
    
    backtrackingOptimizer.apply (problems, objectiveFunctional, g, updater);
            

    // compute difference between target and reconstructed solution
    dolfin::Function difference (V);
    dolfin::Function target (V);
    target = u_0;
    difference = problems.solution ("primal") - target; 
    
    
    // ===============
    // POST PROCESSING
    // ===============
    dolfin::VTKPlotter meshPlotter (dolfin::reference_to_no_delete_pointer (mesh));
    meshPlotter.parameters["title"] = "Mesh";
    meshPlotter.plot ();
    
    dolfin::VTKPlotter solutionPlotter (dolfin::reference_to_no_delete_pointer (problems.solution ("primal")));
    solutionPlotter.parameters["title"] = "Solution of the problem with computed control";
    solutionPlotter.parameters["input_keys"] = "m";
    solutionPlotter.plot ();
    
    dolfin::VTKPlotter controlPlotter (dolfin::reference_to_no_delete_pointer (g));
    controlPlotter.parameters["title"] = "Control";
    controlPlotter.plot ();
    
    dolfin::VTKPlotter controlRegionPlotter (dolfin::reference_to_no_delete_pointer (meshCells));
    controlRegionPlotter.parameters["mode"]  = "color";
    controlRegionPlotter.parameters["title"] = "Control region";
    controlRegionPlotter.plot ();
    
    dolfin::VTKPlotter differencePlotter (dolfin::reference_to_no_delete_pointer (difference));
    differencePlotter.parameters["title"] = "Difference between target and reconstructed solution";
    differencePlotter.parameters["mode"]  = "color";
    differencePlotter.plot ();
    
    dolfin::interactive ();
    
    return 0;
}
