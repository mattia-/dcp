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
    std::cout << "Create mesh and finite element space..." << std::endl;
    dolfin::UnitSquareMesh mesh (100, 100);
    poisson::FunctionSpace V (mesh);
    
    // define problem
    std::cout << "Define the problem..." << std::endl;
    dcp::LinearDifferentialProblem <poisson::BilinearForm, poisson::LinearForm> 
        poissonProblem (dolfin::reference_to_no_delete_pointer (mesh), 
                        dolfin::reference_to_no_delete_pointer (V));

    // define coefficients
    std::cout << "Define the problem's coefficients..." << std::endl;
    dolfin::Constant k (1.0);
    dolfin::Constant f (1.0);
    dolfin::Constant g (1.0);
    dolfin::Constant h (0.0);

    // define dirichlet boundary conditions 
    std::cout << "Define the problem's Dirichlet boundary conditions..." << std::endl;
    poisson::DirichletBoundary dirichletBoundary;
    poissonProblem.addDirichletBC (dolfin::DirichletBC (V, h, dirichletBoundary));
    
    // define neumann boundary conditions 
    std::cout << "Define the problem's Neumann boundary conditions..." << std::endl;
    poisson::NeumannBoundary neumannBoundary;
    dolfin::FacetFunction<std::size_t> meshFacets (mesh);
    meshFacets.set_all (0);
    neumannBoundary.mark (meshFacets, 1);
    
    // set problem coefficients
    std::cout << "Set the problem's coefficients..." << std::endl;
    poissonProblem.setCoefficient ("bilinear_form", dolfin::reference_to_no_delete_pointer (k), "k");
    poissonProblem.setCoefficient ("linear_form", dolfin::reference_to_no_delete_pointer (k), "k");
    poissonProblem.setCoefficient ("linear_form", dolfin::reference_to_no_delete_pointer (f), "f");
    poissonProblem.setCoefficient ("linear_form", dolfin::reference_to_no_delete_pointer (g), "g");

    // solve problem
    std::cout << "Solve the problem..." << std::endl;
    poissonProblem.solve ();
    
    // plots
    dolfin::plot (mesh);
    dolfin::plot (poissonProblem.solution ());
    
    dolfin::interactive ();
    
    return 0;
}
