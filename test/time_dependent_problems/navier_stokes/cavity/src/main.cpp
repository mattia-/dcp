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
#include "navierstokes.h"
#include <dcp/differential_problems/differential_problems.h>

namespace navierstokes
{
    class MovingLid : public dolfin::SubDomain
    {
        bool inside (const dolfin::Array<double>& x, bool on_boundary) const
        {
            return dolfin::near (x[1], 1) && on_boundary;
        }
    };

    class FixedWalls : public dolfin::SubDomain
    {
        bool inside (const dolfin::Array<double>& x, bool on_boundary) const
        {
            return (dolfin::near (x[1], 0) && on_boundary)
                   ||
                   (dolfin::near (x[0], 0) && on_boundary) 
                   ||
                   (dolfin::near (x[0], 1) && on_boundary);
                   
        }
    };
}

int main (int argc, char* argv[])
{
    // create mesh and finite element space
    std::cout << "Create mesh and finite element space..." << std::endl;
    dolfin::UnitSquareMesh mesh (20, 20);
    
    navierstokes::FunctionSpace V (mesh);
    
    // define problem
    std::cout << "Define the problem..." << std::endl;
    
    double t0 = 0.0;
    double dt = 0.1;
    double T = 4;
    dcp::NonlinearProblem <navierstokes::ResidualForm, navierstokes::JacobianForm> 
        timeSteppingProblem (dolfin::reference_to_no_delete_pointer (V), "trial");
    
    dcp::TimeDependentProblem navierStokesProblem (dolfin::reference_to_no_delete_pointer (timeSteppingProblem),
                                                   t0,
                                                   dt, 
                                                   T, 
                                                   {"residual_form", "jacobian_form"},
                                                   {"residual_form"});
         
    // define constant
    std::cout << "Define the problem's coefficients..." << std::endl;
    dolfin::Constant nu (1e-1);
    dolfin::Constant movingLidVelocity (1.0, 0.0);
    dolfin::Constant noSlipDirichletBC (0.0, 0.0);

    // define dirichlet boundary conditions
    std::cout << "Define the problem's Dirichlet boundary conditions..." << std::endl;
    navierstokes::MovingLid movingLid;
    navierstokes::FixedWalls fixedWalls;
    
    navierStokesProblem.addDirichletBC (movingLidVelocity, movingLid, 0);
    navierStokesProblem.addDirichletBC (noSlipDirichletBC, fixedWalls, 0);

    // problem settings
    std::cout << "Set the problem's coefficients..." << std::endl;
    navierStokesProblem.setCoefficient ("residual_form", dolfin::reference_to_no_delete_pointer (nu), "nu");
    navierStokesProblem.setCoefficient ("jacobian_form", dolfin::reference_to_no_delete_pointer (nu), "nu");
    
    navierStokesProblem.parameters ["time_stepping_solution_component"] = 0;
    navierStokesProblem.parameters ["plot_component"] = 0;
    navierStokesProblem.parameters ["plot_interval"] = 1;
    
    // plot mesh
    dolfin::plot (mesh, "Mesh");
    
    // solve problem
    std::cout << "Solve the problem..." << std::endl;
    navierStokesProblem.solve ();

    // dolfin::interactive ();
    
    return 0;
}
