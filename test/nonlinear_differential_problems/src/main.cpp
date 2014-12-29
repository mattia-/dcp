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
#include "navierstokes.h"
#include <DifferentialProblem/NonlinearDifferentialProblem.hpp>

namespace navierstokes
{
    class InflowBoundary : public dolfin::SubDomain
    {
        bool inside (const dolfin::Array<double>& x, bool on_boundary) const
        {
            return x[0] <= (0 + DOLFIN_EPS) && on_boundary;
        }
    };

    class GammaSD : public dolfin::SubDomain
    {
        bool inside (const dolfin::Array<double>& x, bool on_boundary) const
        {
            return (x[1] <= (0 + DOLFIN_EPS) || x[1] >= (7 - DOLFIN_EPS)) && on_boundary;
        }
    };
    
    class Circle : public dolfin::SubDomain
    {
        public:
            Circle () : 
                center (2),
                radius (0.5)
            {
                center[0] = 3.5;
                center[1] = 3.5;
            }
            
            bool inside (const dolfin::Array<double>& x, bool on_boundary) const
            {
                double dx = x[0] - center[0];
                double dy = x[1] - center[1];
                double r = sqrt (dx * dx + dy * dy);
                
                return r <= (radius + 1e-3) && on_boundary;
            }

        private:
            dolfin::Array<double> center;
            double radius;
    }; 
    
    class NoSlipBoundary : public dolfin::SubDomain
    { 
        public:
            bool inside (const dolfin::Array<double>& x, bool on_boundary) const
            {
                bool onWall1 = (x[0] >= 1.75 - DOLFIN_EPS && x[0] <= 1.75 + DOLFIN_EPS)
                            && (x[1] >= 3 - DOLFIN_EPS && x[1] <= 4 + DOLFIN_EPS)
                            && on_boundary;
                bool onWall2 = (x[0] >= (1.75 - DOLFIN_EPS) && x[0] <= (2 + DOLFIN_EPS))
                            && (x[1] <= (4 + DOLFIN_EPS) && x[1] >= (4 - DOLFIN_EPS))
                            && on_boundary;
                bool onWall3 = (x[0] >= 2 - DOLFIN_EPS && x[0] <= 2 + DOLFIN_EPS)
                            && (x[1] >= 3 - DOLFIN_EPS && x[1] <= 4 + DOLFIN_EPS)
                            && on_boundary;
                bool onWall4 = (x[0] >= (1.75 - DOLFIN_EPS) && x[0] <= (2 + DOLFIN_EPS))
                            && (x[1] <= (3 + DOLFIN_EPS) && x[1] >= (3 - DOLFIN_EPS))
                            && on_boundary;
                
                bool onCircle = circle.inside (x, on_boundary);

                return onWall1 || onWall2 || onWall3 || onWall4 || onCircle;
            }

        private:
            navierstokes::Circle circle;
    }; 
}

int main (int argc, char* argv[])
{
    // create mesh and finite element space
    std::cout << "Create mesh and finite element space..." << std::endl;
    dolfin::Mesh mesh ("../../optimization/src/complete_mesh/complete_mesh.xml");
    navierstokes::FunctionSpace V (mesh);
    
    // define problem
    std::cout << "Define the problem..." << std::endl;
    dcp::NonlinearDifferentialProblem <navierstokes::ResidualForm, navierstokes::JacobianForm> 
        navierStokesProblem (dolfin::reference_to_no_delete_pointer (mesh), 
                             dolfin::reference_to_no_delete_pointer (V),
                             "trial");
    
    // define constant
    std::cout << "Define the problem's coefficients..." << std::endl;
    dolfin::Constant nu (1e-1);
    dolfin::Constant inflowDirichletBC (1.0, 0.0);
    dolfin::Constant symmetryDirichletBC (0.0);
    dolfin::Constant noSlipDirichletBC (0.0, 0.0);

    // define dirichlet boundary conditions
    std::cout << "Define the problem's Dirichlet boundary conditions..." << std::endl;
    navierstokes::InflowBoundary inflowBoundary;
    navierstokes::GammaSD gammaSD;
    navierstokes::NoSlipBoundary noSlipBoundary;
    
    navierStokesProblem.addDirichletBC (dolfin::DirichletBC (*V[0], inflowDirichletBC, inflowBoundary));
    navierStokesProblem.addDirichletBC (dolfin::DirichletBC (*V[0], noSlipDirichletBC, noSlipBoundary, "topological", false));
    navierStokesProblem.addDirichletBC (dolfin::DirichletBC (*(*V[0])[1], symmetryDirichletBC, gammaSD));

    // problem settings
    std::cout << "Set the problem's coefficients..." << std::endl;
    navierStokesProblem.setCoefficient ("residual_form", dolfin::reference_to_no_delete_pointer (nu), "nu");
    navierStokesProblem.setCoefficient ("jacobian_form", dolfin::reference_to_no_delete_pointer (nu), "nu");
    
    // solve problem
    std::cout << "Solve the problem..." << std::endl;
    navierStokesProblem.solve ();

    // plots
    dolfin::plot (mesh, "Mesh");
    dolfin::plot (navierStokesProblem.solution ()[0], "Velocity");
    dolfin::plot (navierStokesProblem.solution ()[1], "Pressure");
    dolfin::interactive ();
    
    return 0;
}
