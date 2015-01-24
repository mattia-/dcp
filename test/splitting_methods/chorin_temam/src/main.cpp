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
#include "projection.h"
#include <differential_problems/differential_problems.h>
#include <splitting_methods/splitting_methods.h>

namespace navierstokes
{
    class InflowBoundary : public dolfin::SubDomain
    {
        bool inside (const dolfin::Array<double>& x, bool on_boundary) const
        {
            return x[0] <= (0 + DOLFIN_EPS) && on_boundary;
        }
    };

    class OutflowBoundary : public dolfin::SubDomain
    {
        bool inside (const dolfin::Array<double>& x, bool on_boundary) const
        {
            return x[0] >= (10 - DOLFIN_EPS) && on_boundary;
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
    mshr::Rectangle rectangle (dolfin::Point (0.0, 0.0), dolfin::Point (10.0, 7.0));
    mshr::Circle circle (dolfin::Point (3.5, 3.5), 0.5);
    dolfin::Mesh mesh;
    mshr::generate (mesh, *(rectangle - circle), 50);
    
    prediction::FunctionSpace V (mesh);
    correction::FunctionSpace Q (mesh);
    
    // define constant
    std::cout << "Define the coefficients..." << std::endl;
    dolfin::Constant nu (1e-1);
    dolfin::Constant inflowDirichletBC (1.0, 0.0);
    dolfin::Constant symmetryDirichletBC (0.0);
    dolfin::Constant noSlipDirichletBC (0.0, 0.0);
    double t0 = 0.0;
    double dt = 0.1;
    double T = 4;
    
    // define the problems
    std::cout << "Define the problem..." << std::endl;
    dcp::ChorinTemamMethod <prediction::BilinearForm,
                            prediction::LinearForm,
                            correction::BilinearForm,
                            correction::LinearForm,
                            projection::BilinearForm,
                            projection::LinearForm>
        chorinTemamMethod ({V, Q}, t0, dt, T, nu);                                
    
    
    // define dirichlet boundary conditions
    std::cout << "Define the problem's Dirichlet boundary conditions..." << std::endl;
    navierstokes::InflowBoundary inflowBoundary;
    navierstokes::OutflowBoundary outflowBoundary;
    navierstokes::GammaSD gammaSD;
    navierstokes::NoSlipBoundary noSlipBoundary;
    
    std::cout << "Setting up dirichlet's boundary conditions for 'prediction_problem'" << std::endl;
    chorinTemamMethod.addDirichletBC ("prediction_problem", inflowDirichletBC, inflowBoundary);
    chorinTemamMethod.addDirichletBC ("prediction_problem", noSlipDirichletBC, noSlipBoundary);
    chorinTemamMethod.addDirichletBC ("prediction_problem", symmetryDirichletBC, gammaSD, 1);

    std::cout << "Setting up dirichlet's boundary conditions for 'correction_problem'" << std::endl;
    chorinTemamMethod.addDirichletBC ("correction_problem", symmetryDirichletBC, gammaSD);
    chorinTemamMethod.addDirichletBC ("correction_problem", symmetryDirichletBC, outflowBoundary);
    
    for (auto& i : chorinTemamMethod.system ().problemsNames ())
    {
        chorinTemamMethod [i].parameters ["plot_interval"] = 1;
        chorinTemamMethod [i].parameters ["plot_title"] = i;
    }
    
    // plots
    dolfin::plot (mesh, "Mesh");
    
    // solve problem
    std::cout << "Solve the problem..." << std::endl;
    chorinTemamMethod.apply ();

    dolfin::interactive ();
    
    return 0;
}
