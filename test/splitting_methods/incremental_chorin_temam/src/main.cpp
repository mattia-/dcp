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
#include "prediction.h"
#include "correction.h"
#include "velocity_projection.h"
#include "pressure_update.h"
#include <dcp/differential_problems/differential_problems.h>
#include <dcp/splitting_methods/splitting_methods.h>

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

    prediction::FunctionSpace V (mesh);
    correction::FunctionSpace Q (mesh);

    // define constant
    std::cout << "Define the coefficients..." << std::endl;
    double t0 = 0.0;
    double dt = 0.1;
    double T = 4;
    dolfin::Constant nu (1e-1);
    dolfin::Constant movingLidVelocity (1.0, 0.0);
    dolfin::Constant noSlipDirichletBC (0.0, 0.0);

    // define the problems
    std::cout << "Define the problem..." << std::endl;
    dcp::IncrementalChorinTemamMethod <prediction::BilinearForm,
                                       prediction::LinearForm,
                                       correction::BilinearForm,
                                       correction::LinearForm,
                                       velocity_projection::BilinearForm,
                                       velocity_projection::LinearForm,
                                       pressure_update::BilinearForm,
                                       pressure_update::LinearForm>
        incrementalChorinTemamMethod ({V, Q}, t0, dt, T, nu);


    // define dirichlet boundary conditions
    std::cout << "Define the problem's Dirichlet boundary conditions..." << std::endl;
    navierstokes::MovingLid movingLid;
    navierstokes::FixedWalls fixedWalls;

    std::cout << "Setting up dirichlet's boundary conditions for \"prediction_problem\"" << std::endl;
    incrementalChorinTemamMethod.addDirichletBC ("prediction_problem", movingLidVelocity, movingLid);
    incrementalChorinTemamMethod.addDirichletBC ("prediction_problem", noSlipDirichletBC, fixedWalls);

    for (auto& i : incrementalChorinTemamMethod.system ().problemsNames ())
    {
        incrementalChorinTemamMethod [i].parameters ["plot_interval"] = 1;
        incrementalChorinTemamMethod [i].parameters ["plot_title"] = i;
    }
    
    // plots
    dolfin::plot (mesh, "Mesh");

    // solve problem
    std::cout << "Solve the problem..." << std::endl;
    incrementalChorinTemamMethod.apply ();

    dolfin::interactive ();

    return 0;
}
