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
#include <dcp/expressions/expressions.h>
#include <dcp/subdomains/subdomains.h>

namespace navierstokes
{
    // the external force is divided into two parts: f = f1 * rho + f2
    // this is because rho is a solution, and there is no way to link it
    // to a TimeDependentExpression
    class FirstExternalLoad : public dcp::TimeDependentExpression
    {
        public:
            FirstExternalLoad () : TimeDependentExpression (2) { }
        
            void eval (dolfin::Array<double>& values, const dolfin::Array<double>& x) const override
            {
                values[0] =  x[1] * sin (t) - x[0] * cos (t) * cos (t);
                values[1] = -x[0] * sin (t) - x[1] * cos (t) * cos (t);
            }
    }; 
    

    class SecondExternalLoad : public dcp::TimeDependentExpression
    {
        public:
            SecondExternalLoad () : TimeDependentExpression (2) { }

            void eval (dolfin::Array<double>& values, const dolfin::Array<double>& x) const override
            {
                values[0] =  cos (x[0]) * sin (x[1]) * sin (t);
                values[1] =  sin (x[0]) * cos (x[1]) * sin (t);
            }
    }; 
    
    
    class DirichletBoundaryEvaluator
    { 
        public:
            bool operator() (const dolfin::Array<double>& x, bool on_boundary)
            {
                return on_boundary;
            }
    }; 
}

namespace exact_solutions
{
    class RhoEvaluator
    {
        public:
            void operator() (dolfin::Array<double>& values, const dolfin::Array<double>& x, const double& t)
            {
                values[0] = 2 + x[0] * cos (sin (t)) + x[1] * sin (sin (t));
            }
    };
    
    class UEvaluator
    {
        public:
            void operator() (dolfin::Array<double>& values, const dolfin::Array<double>& x, const double& t)
            {
                values[0] = -x[1] * cos (t);
                values[1] =  x[0] * cos (t);
            }
    };
    
    class PEvaluator
    {
        public:
            void operator() (dolfin::Array<double>& values, const dolfin::Array<double>& x, const double& t)
            {
                values[0] = sin (x[0]) * sin (x[1]) * sin (t);
            }
    };
}

int main (int argc, char* argv[])
{
//    dolfin::set_log_level (dolfin::DBG);
    // create mesh and finite element space
    std::cout << "Create mesh and finite element space..." << std::endl;
    mshr::Circle circle (dolfin::Point (0.0, 0.0), 1);
    dolfin::Mesh mesh;
    mshr::generate (mesh, circle, 50);
    
    // define constant
    std::cout << "Define the coefficients..." << std::endl;
    dolfin::Constant mu (1e-1);
    dolfin::Constant chi (1);
    double t0 = 0.0;
    double dt = 0.1;
    double T = 10;
    
    // exact solution
    dcp::TimeDependentExpression exactRho ((exact_solutions::RhoEvaluator ()));
    dcp::TimeDependentExpression exactU (2, exact_solutions::UEvaluator ());
    dcp::TimeDependentExpression exactP ((exact_solutions::PEvaluator ()));

    double t = t0;
    while (t < T)
    {
        exactRho.setTime (t);
        exactU.setTime (t);
        exactP.setTime (t);
        
        dolfin::plot (exactRho, mesh, "rho: time = " + std::to_string (t));
        dolfin::plot (exactU, mesh, "u: time = " + std::to_string (t));
        dolfin::plot (exactP, mesh, "p: time = " + std::to_string (t));

        dolfin::interactive ();
        
        t += dt;
    }
    
    return 0;
}

