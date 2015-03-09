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
#include <differential_problems/differential_problems.h>
#include "laplaceVec.h"

namespace laplaceVec
{
/*    class NeumannBoundary : public dolfin::SubDomain
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
*/
    class FreeSurface : public dolfin::SubDomain
    {
        bool inside (const dolfin::Array<double>& x, bool on_boundary) const
        {
            return on_boundary && 
                   !( dolfin::near (x[0], 0)
                   ||
                   dolfin::near (x[1], 0)
                   ||
                   dolfin::near (x[0], 1));
        }
/*        bool inside (const dolfin::Array<double>& x, bool on_boundary) const
        {
            return dolfin::near (x[1],1) && on_boundary;
        }*/
    };
    
    class FixedBoundary : public dolfin::SubDomain
    {
        bool inside (const dolfin::Array<double>& x, bool on_boundary) const
        {
            return (dolfin::near (x[0], 0) && on_boundary)
                   ||
                   (dolfin::near (x[1], 0) && on_boundary)
                   ||
                   (dolfin::near (x[0], 1) && on_boundary);
        }
    };
    
    class ParabolicProfile : public dolfin::Expression
    {
        void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x) const
        {
            values[0] = 0.5*x[0]*(1-x[0]);
            values[1] = 0.5*x[0]*(1-x[0]);
        }
        
        std::size_t value_rank() const
        {
          return 1;
        }
    
        std::size_t value_dimension(std::size_t i) const
        {
          return 2;
        }
    };
}

int main (int argc, char* argv[])
{
    // create mesh and finite element space 
    dolfin::info ("Create mesh and finite element space...");
    dolfin::UnitSquareMesh mesh (20, 20);
    laplaceVec::FunctionSpace V (mesh);
    
    // define problem
    dolfin::info ("Define the problem...");
    dcp::LinearProblem <laplaceVec::BilinearForm, laplaceVec::LinearForm> laplaceVecProblem (dolfin::reference_to_no_delete_pointer (mesh), 
                                                                                    dolfin::reference_to_no_delete_pointer (V));

    // define coefficients
    dolfin::info ("Define the problem's coefficients...");
    dolfin::Constant zero (0,0);
    dolfin::Constant uno (1,1);

    // define dirichlet boundary conditions 
    dolfin::info ("Define the problem's Dirichlet boundary conditions...");
    laplaceVec::FreeSurface freeSurface;
    laplaceVec::FixedBoundary fixedBoundary;
    laplaceVec::ParabolicProfile parabolicProfile;
//    laplaceVecProblem.addDirichletBC (dolfin::DirichletBC (V, parabolicProfile, freeSurface));
    laplaceVecProblem.addDirichletBC (dolfin::DirichletBC (V, uno, freeSurface));
    laplaceVecProblem.addDirichletBC (dolfin::DirichletBC (V, zero, fixedBoundary));
    
/*    // define neumann boundary conditions 
    dolfin::info ("Define the problem's Neumann boundary conditions...");
    laplaceVec::NeumannBoundary neumannBoundary;
    dolfin::FacetFunction<std::size_t> meshFacets (mesh);
    meshFacets.set_all (0);
    neumannBoundary.mark (meshFacets, 1);
*/
    
    // set problem coefficients
/*    dolfin::info ("Set the problem's coefficients...");
    laplaceVecProblem.setCoefficient ("bilinear_form", dolfin::reference_to_no_delete_pointer (k), "k");
    laplaceVecProblem.setCoefficient ("linear_form", dolfin::reference_to_no_delete_pointer (k), "k");
    laplaceVecProblem.setCoefficient ("linear_form", dolfin::reference_to_no_delete_pointer (f), "f");
    laplaceVecProblem.setCoefficient ("linear_form", dolfin::reference_to_no_delete_pointer (g), "g");*/
    laplaceVecProblem.setCoefficient ("linear_form", dolfin::reference_to_no_delete_pointer (zero), "zero");

    // solve problem
    dolfin::info ("Solve the problem...");
    laplaceVecProblem.solve ();
    
    // plots
    dolfin::plot (mesh);
    dolfin::plot (laplaceVecProblem.solution ());
    
    dolfin::interactive ();
    
    return 0;
}
