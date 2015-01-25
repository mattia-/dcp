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
#include <dcp/differential_problems/differential_problems.h>
#include "poisson.h"

namespace poisson
{
    class NeumannBoundary : public dolfin::SubDomain
    {
        bool inside (const dolfin::Array<double>& x, bool on_boundary) const
        {
            return dolfin::near (x[1], 0) && on_boundary;
        } 
    };
    
    class DirichletBoundary : public dolfin::SubDomain
    {
        bool inside (const dolfin::Array<double>& x, bool on_boundary) const
        {
            return (dolfin::near (x[0], 0) && on_boundary)
                   ||
                   (dolfin::near (x[1], 1) && on_boundary)
                   ||
                   (dolfin::near (x[0], 1) && on_boundary);
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
    dolfin::Constant g (1.0);
    dolfin::Constant h (0.0);

    // define dirichlet boundary conditions 
    dolfin::info ("Define the problem's Dirichlet boundary conditions...");
    poisson::DirichletBoundary dirichletBoundary;
    poissonProblem.addDirichletBC (h, dirichletBoundary);
    
    // define neumann boundary conditions 
    dolfin::info ("Define the problem's Neumann boundary conditions...");
    poisson::NeumannBoundary neumannBoundary;
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
    
    dolfin::interactive ();
    
    return 0;
}
