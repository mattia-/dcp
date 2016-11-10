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
#include <dcp.h>

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
    dolfin::set_log_level (dolfin::DBG);

    // create mesh and finite element space
    std::cout << "Create mesh and finite element space..." << std::endl;
    mshr::Rectangle rectangle (dolfin::Point (0.0, 0.0), dolfin::Point (10.0, 7.0));
    mshr::Circle circle (dolfin::Point (3.5, 3.5), 0.5);
    auto mesh = mshr::generate_mesh (*(rectangle - circle), 50);
    
    auto V = std::make_shared<navierstokes::FunctionSpace> (mesh);
    
    // define problem
    std::cout << "Define the problem..." << std::endl;
    dcp::NonlinearProblem <navierstokes::ResidualForm, navierstokes::JacobianForm> 
        navierStokesProblem (V, "trial");
    
    // define constant
    std::cout << "Define the problem's coefficients..." << std::endl;
    auto nu = std::make_shared<dolfin::Constant> (1e-1);
    auto inflowDirichletBC = std::make_shared<dolfin::Constant> (1.0, 0.0);
    auto symmetryDirichletBC = std::make_shared<dolfin::Constant> (0.0);
    auto noSlipDirichletBC = std::make_shared<dolfin::Constant> (0.0, 0.0);

    // define dirichlet boundary conditions
    std::cout << "Define the problem's Dirichlet boundary conditions..." << std::endl;
    auto inflowBoundary = std::make_shared<navierstokes::InflowBoundary> ();
    auto gammaSD = std::make_shared<navierstokes::GammaSD> ();
    auto noSlipBoundary = std::make_shared<navierstokes::NoSlipBoundary> ();
    
    navierStokesProblem.addDirichletBC (inflowDirichletBC, inflowBoundary, 0);
    navierStokesProblem.addDirichletBC (dolfin::DirichletBC ((*V)[0], noSlipDirichletBC, noSlipBoundary, "topological", false));
    navierStokesProblem.addDirichletBC (dolfin::DirichletBC ((*(*V)[0])[1], symmetryDirichletBC, gammaSD));

    // problem settings
    std::cout << "Set the problem's coefficients..." << std::endl;
    navierStokesProblem.setCoefficient ("residual_form", nu, "nu");
    navierStokesProblem.setCoefficient ("jacobian_form", nu, "nu");
    
    // solve problem
    std::cout << "Solve the problem..." << std::endl;
    navierStokesProblem.solve ();

    // plots
    dolfin::plot (mesh, "Mesh");
    navierStokesProblem.parameters ["plot_components"] = "0 1";
    navierStokesProblem.plotSolution ();
    /* ALSO:
     * dolfin::plot (navierStokesProblem.solution ()[0], "Velocity");
     * dolfin::plot (navierStokesProblem.solution ()[1], "Pressure");
     */
    // dolfin::interactive ();
    
    return 0;
}
