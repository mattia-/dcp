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
#include "convectiondiffusion.h"

namespace convectiondiffusion
{
    class DirichletBoundary : public dolfin::SubDomain
    {
        bool inside (const dolfin::Array<double>& x, bool on_boundary) const
        {
            return on_boundary;
        }
    };
    
    class initialSolution : public dolfin::Expression
    {
        void eval (dolfin::Array<double>& values, const dolfin::Array<double>& x) const
        {
            values[0] = exp (-((x[0] - 1) * (x[0] - 1) + (x[1] - 3) * (x[1] - 3)) / 0.2);
        }
    };
}

int main (int argc, char* argv[])
{
    // create mesh and finite element space 
    std::cout << "Create mesh and finite element space..." << std::endl;
    mshr::Rectangle rectangle (dolfin::Point (0.0, 0.0), dolfin::Point (10.0, 7.0));
    mshr::Circle circle (dolfin::Point (3.5, 3.5), 0.5);
    dolfin::Mesh mesh;
    mshr::generate (mesh, *(rectangle - circle), 50);
    
    convectiondiffusion::FunctionSpace V (mesh);
    
    // define problem
    std::cout << "Define the problem..." << std::endl;
    double t0 = 0;
    double dt = 0.1;
    double T = 4;
    
    dcp::LinearProblem <convectiondiffusion::BilinearForm, convectiondiffusion::LinearForm>
        timeSteppingProblem (dolfin::reference_to_no_delete_pointer (V));
    
    dcp::TimeDependentProblem convectionDiffusionProblem (dolfin::reference_to_no_delete_pointer (timeSteppingProblem),
                                                          t0,
                                                          dt,
                                                          T,
                                                          {"bilinear_form"},
                                                          {"linear_form"});
    
    
    // define coefficients
    std::cout << "Define the problem's coefficients..." << std::endl;
    dolfin::Constant k (0.01);
    dolfin::Constant b (1.0, 0.25); 

    // define dirichlet boundary conditions 
    std::cout << "Define the problem's Dirichlet boundary conditions..." << std::endl;
    convectiondiffusion::DirichletBoundary dirichletBoundary;
    dolfin::Constant dirichletValue (0.0);
    convectionDiffusionProblem.addDirichletBC (dirichletValue, dirichletBoundary);
    
    // set problem coefficients
    std::cout << "Set the problem's coefficients..." << std::endl;
    convectionDiffusionProblem.setCoefficient ("bilinear_form", dolfin::reference_to_no_delete_pointer (k), "k");
    convectionDiffusionProblem.setCoefficient ("bilinear_form", dolfin::reference_to_no_delete_pointer (b), "b");

    // solve problem
    convectionDiffusionProblem.setInitialSolution (convectiondiffusion::initialSolution ());
    convectionDiffusionProblem.parameters ["plot_interval"] = 1;
    std::cout << "Solve the problem..." << std::endl;
    convectionDiffusionProblem.solve ();
  
    // plots
    dolfin::plot (mesh);
    
    dolfin::interactive ();
  
    return 0;
}
