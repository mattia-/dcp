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
#include <dcp/problems/problems.h>
#include <dcp/expressions/expressions.h>
#include <dcp/subdomains/subdomains.h>
#include "poisson.h"

namespace poisson
{
    class NeumannBoundary
    {
        public:
            bool operator() (const dolfin::Array<double>& x, bool on_boundary)
            {
                return dolfin::near (x[1], 0) && on_boundary;
            } 
    };
    
    class DirichletBoundary
    {
        public:
            bool operator() (const dolfin::Array<double>& x, bool on_boundary)
            {
                return (dolfin::near (x[0], 0) && on_boundary)
                       ||
                       (dolfin::near (x[1], 1) && on_boundary)
                       ||
                       (dolfin::near (x[0], 1) && on_boundary);
            }
    };
    
    class NeumannExpression
    {
        public:
            void operator() (dolfin::Array<double>& values, const dolfin::Array<double>& x)
            {
                values[0] = x[0] * x[0];
            }
    };
    
    class DirichletExpression
    {
        public:
            void operator() (dolfin::Array<double>& values, const dolfin::Array<double>& x)
            {
                values[0] = x[0] + x[1];
            }
    };
}

int main (int argc, char* argv[])
{
    // create mesh and finite element space 
    dolfin::info ("Create mesh and finite element space...");
    dolfin::UnitSquareMesh mesh (20, 20);
    poisson::FunctionSpace V (mesh);
    
    // define problem
    dolfin::info ("Define the problem...");
    dcp::LinearProblem <poisson::BilinearForm, poisson::LinearForm> poissonProblem (dolfin::reference_to_no_delete_pointer (V));

    // define coefficients
    dolfin::info ("Define the problem's coefficients...");
    dolfin::Constant k (1.0);
    dolfin::Constant f (1.0);
    dcp::Expression g ((poisson::NeumannExpression ()));
    dcp::Expression h ((poisson::DirichletExpression ()));

    // define dirichlet boundary conditions 
    dolfin::info ("Define the problem's Dirichlet boundary conditions...");
    dcp::Subdomain dirichletBoundary ((poisson::DirichletBoundary ()));
    poissonProblem.addDirichletBC (h, dirichletBoundary);
    
    // define neumann boundary conditions 
    dolfin::info ("Define the problem's Neumann boundary conditions...");
    dcp::Subdomain neumannBoundary ((poisson::NeumannBoundary ()));
    dolfin::FacetFunction<std::size_t> meshFacets (mesh);
    meshFacets.set_all (0);
    neumannBoundary.mark (meshFacets, 1);
    
    // set problem coefficients
    dolfin::info ("Set the problem's coefficients...");
    poissonProblem.setCoefficient ("bilinear_form", dolfin::reference_to_no_delete_pointer (k), "k");
    poissonProblem.setCoefficient ("linear_form", dolfin::reference_to_no_delete_pointer (k), "k");
    poissonProblem.setCoefficient ("linear_form", dolfin::reference_to_no_delete_pointer (f), "f");
    poissonProblem.setCoefficient ("linear_form", dolfin::reference_to_no_delete_pointer (g), "g");

    // solve problem
    dolfin::info ("Solve the problem...");
    poissonProblem.solve ();
    
    // plots
    dolfin::plot (mesh);
    poissonProblem.plotSolution ();
    
    // dolfin::interactive ();
    
    return 0;
}
