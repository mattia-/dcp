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
#include <dcp.h>
#include "primal.h"
#include "adjoint.h"
#include "objective_functional.h"

class DirichletBoundary
{
    public:
        bool operator() (const dolfin::Array<double>& x, bool on_boundary)
        {
            return on_boundary;
        }
};

class TargetSolutionEvaluator
{
    public:
        void operator() (dolfin::Array<double>& values, const dolfin::Array<double>& x)
        {
            values[0] = x[0] * (1 - x[0]) * x[1] * (1 - x[1]);
        }
};

class GradientEvaluator
{
    public:
        void operator() (dolfin::Array<double>& values,
                         const dolfin::Array<double>& x,
                         const std::map <std::string, std::shared_ptr<const dolfin::GenericFunction> >& variables)
        {
            variables.find("u_hat")->second -> eval (values, x);
        }
};


int main (int argc, char* argv[])
{
    dolfin::set_log_level (dolfin::PROGRESS);

    // mesh
    auto mesh = std::make_shared<dolfin::UnitSquareMesh> (20, 20);

    // function spaces
    auto V = std::make_shared<primal::FunctionSpace> (mesh);
    auto W = std::make_shared<adjoint::FunctionSpace> (mesh);


    // =======================
    // PRIMAL+ADJOINT PROBLEMS
    // =======================
    // define problems
    dcp::LinearProblem <primal::BilinearForm, primal::LinearForm> primalProblem (V);

    dcp::LinearProblem <adjoint::BilinearForm, adjoint::LinearForm> adjointProblem (W);

    // create composite differential problem
    dcp::EquationSystem problems;
    problems.addProblem ("primal", primalProblem);
    problems.addProblem ("adjoint", adjointProblem);

    // define constants
    dcp::Expression u_0 ((TargetSolutionEvaluator ()));
    dolfin::Constant dirichletBC (0.0);

    // define boundary conditions subdomains
    dcp::Subdomain dirichletBoundary ((DirichletBoundary ()));

    // primal problem settings
    problems["primal"].addDirichletBC (dirichletBC, dirichletBoundary);

    // adjoint problem settings
    problems["adjoint"].addDirichletBC (dirichletBC, dirichletBoundary);
    problems["adjoint"].setCoefficient ("linear_form", dolfin::reference_to_no_delete_pointer (u_0), "u_0");

    problems.addLink ("adjoint", "u", "linear_form", "primal");


    // ====================
    // OBJECTIVE FUNCTIONAL
    // ====================
    // define functional
    dcp::ObjectiveFunctional <objective_functional::Form_J>
        objectiveFunctional (mesh, std::make_shared<dcp::VariableExpression> (GradientEvaluator ()));

    // control variable
    dolfin::Function g (V);
    g = dolfin::Constant (0.0);

    // functional settings
    objectiveFunctional.setCoefficient ("functional", dolfin::reference_to_no_delete_pointer (u_0), "u_0");
    objectiveFunctional.setCoefficient ("functional",
                                        dolfin::reference_to_no_delete_pointer (problems.solution ("primal")),
                                        "u");

    objectiveFunctional.setCoefficient ("gradient",
                                        dolfin::reference_to_no_delete_pointer (problems.solution ("adjoint")),
                                        "u_hat");


    // ============
    // OPTIMIZATION
    // ============
    problems["primal"].setCoefficient ("linear_form", dolfin::reference_to_no_delete_pointer (g), "g");

    // define optimizer
    dcp::BacktrackingOptimizer backtrackingOptimizer;
    backtrackingOptimizer.parameters ["output_file_name"] = "results.txt";
    backtrackingOptimizer.parameters ["max_minimization_iterations"] = 100;

    dcp::BacktrackingImplementer<dolfin::Function> implementer
        (dcp::DistributedControlUpdater ("primal", "linear_form", "g"));

    backtrackingOptimizer.apply (problems, objectiveFunctional, g, implementer);


    // compute difference between target and reconstructed solution
    dolfin::Function difference (V);
    dolfin::Function target (V);
    target = u_0;
    difference = problems.solution ("primal") - target;


    // ===============
    // POST PROCESSING
    // ===============
    dolfin::plot (mesh, "Mesh");
    dolfin::plot (problems.solution ("primal"), "Primal Solution");
    dolfin::plot (target, "Target Solution");
    dolfin::plot (g, "Control");
    dolfin::plot (difference, "Difference between target and recovered");
    dolfin::plot (mesh, "Mesh");
    dolfin::plot (problems.solution ("primal"), "Primal Solution");
    dolfin::plot (target, "Target Solution");
    dolfin::plot (g, "Control");
    dolfin::plot (difference, "Difference between target and recovered");
    // dolfin::interactive ();

    return 0;
}
