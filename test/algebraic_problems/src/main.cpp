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

class Solver
{
    public:
        void operator() (dolfin::Array<double>& values, 
                         const dolfin::Array<double>& x,
                         const std::map <std::string, std::shared_ptr<const dolfin::GenericFunction> >& variables)
        {
            dolfin::Array<double> uValues (1);
            (*variables.find("u")).second->eval (uValues, x);
            
            if (uValues[0] > 0.25)
            {
                values[0] = 1;
            }
            else
            {
                values[0] = 0;
            }
        }
};

int main (int argc, char* argv[])
{
    // create mesh and finite element space 
    dolfin::info ("Create mesh and finite element space...");
    auto mesh = std::make_shared<dolfin::UnitSquareMesh> (20, 20);
    auto V = std::make_shared<poisson::FunctionSpace> (mesh);
    
    // define linear problem
    dolfin::info ("Define the linear problem...");
    dcp::LinearProblem <poisson::BilinearForm, poisson::LinearForm> poissonProblem (V);

    // define coefficients
    dolfin::info ("Define the linear problem's coefficients...");
    auto k = std::make_shared<dolfin::Constant> (1.0);
    auto f = std::make_shared<dolfin::Constant> (1.0);
    auto g = std::make_shared<dolfin::Constant> (1.0);
    dolfin::Constant h (0.0);

    // define dirichlet boundary conditions 
    dolfin::info ("Define the linear problem's Dirichlet boundary conditions...");
    poisson::DirichletBoundary dirichletBoundary;
    poissonProblem.addDirichletBC (h, dirichletBoundary);
    
    // define neumann boundary conditions 
    dolfin::info ("Define the linear problem's Neumann boundary conditions...");
    poisson::NeumannBoundary neumannBoundary;
    dolfin::FacetFunction<std::size_t> meshFacets (mesh);
    meshFacets.set_all (0);
    neumannBoundary.mark (meshFacets, 1);
    poissonProblem.setIntegrationSubdomain ("linear_form", 
                                             dolfin::reference_to_no_delete_pointer (meshFacets), 
                                             dcp::SubdomainType::BOUNDARY_FACETS);
    
    // set problem coefficients
    dolfin::info ("Set the linear problem's coefficients...");
    poissonProblem.setCoefficient ("bilinear_form", k, "k");
    poissonProblem.setCoefficient ("linear_form", f, "f");
    poissonProblem.setCoefficient ("linear_form", g, "g");

    // solve problem
    dolfin::info ("Solve the linear problem...");
    poissonProblem.solve ();
    
    // plots
    dolfin::plot (mesh);
    poissonProblem.plotSolution ();
    
    // define algebraic problem
    dolfin::info ("Define the algebraic problem...");
    dcp::AlgebraicProblem problem (V, 
                                   std::shared_ptr<dcp::GenericExpression> (new dcp::VariableExpression (Solver ())));
    
    // set algebraic problem coefficients
    dolfin::info ("Set the algebraic problem's coefficients...");
    problem.setCoefficient ("expression", dolfin::reference_to_no_delete_pointer (poissonProblem.solution ()), "u");

    // solve problem
    dolfin::info ("Solve the algebraic problem...");
    problem.solve ();
    
    // plots
    dolfin::plot (mesh);
    problem.plotSolution ();
    
    // dolfin::interactive ();
    
    return 0;
}
