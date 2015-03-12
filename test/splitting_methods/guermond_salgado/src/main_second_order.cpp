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
#include <dcp/differential_problems/differential_problems.h>
#include <dcp/splitting_methods/splitting_methods.h>
#include <dcp/expressions/expressions.h>
#include <dcp/subdomains/subdomains.h>
#include "density_second_order.h"
#include "velocity_second_order.h"
#include "pressure_correction_second_order.h"
#include "pressure_update_second_order.h"
#include "error_computers.h"
#include <fstream>

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
    // get dt and T values from command line. Default value: dt = 0.1, T = 1.0, N = 50
    dolfin::Parameters parameters ("main_parameters");
    parameters.add ("dt", 0.1);
    parameters.add ("T", 1.0);
    parameters.add ("N", 50);
    parameters.parse (argc, argv);
        
//    dolfin::set_log_level (dolfin::DBG);
    // create mesh and finite element space
    std::cout << "Create mesh and finite element space..." << std::endl;
    mshr::Circle circle (dolfin::Point (0.0, 0.0), 1);
    dolfin::Mesh mesh;
    int N = parameters["N"];
    mshr::generate (mesh, circle, N);
    
    density_second_order::FunctionSpace R (mesh);
    velocity_second_order::FunctionSpace V (mesh);
    pressure_correction_second_order::FunctionSpace Q (mesh);
    
    // define constant
    std::cout << "Define the coefficients..." << std::endl;
    dolfin::Constant mu (1e-1);
    dolfin::Constant chi (1);
    double t0 = 0.0;
    double dt = parameters["dt"];
    double T = parameters["T"];
    std::cout << "Time step = " << dt << std::endl; 
    std::cout << "Final time = " << T << std::endl; 
    
    // define the problems
    std::cout << "Define the problem..." << std::endl;
    dcp::GuermondSalgadoMethod <density_second_order::BilinearForm,
                                density_second_order::LinearForm,
                                velocity_second_order::BilinearForm,
                                velocity_second_order::LinearForm,
                                pressure_correction_second_order::BilinearForm,
                                pressure_correction_second_order::LinearForm,
                                pressure_update_second_order::BilinearForm,
                                pressure_update_second_order::LinearForm>
        guermondSalgadoMethod (R, V, Q, t0, dt, T, mu, chi, 2);
    
    // exact solution
    dcp::TimeDependentExpression exactRho ((exact_solutions::RhoEvaluator ()));
    dcp::TimeDependentExpression exactU (2, exact_solutions::UEvaluator ());
    dcp::TimeDependentExpression exactP ((exact_solutions::PEvaluator ()));

    // set initial solutions
    // first step
    exactRho.setTime (t0);
    exactU.setTime (t0);
    exactP.setTime (t0);
    guermondSalgadoMethod.setInitialSolution ("density_problem", exactRho, 2);
    guermondSalgadoMethod.setInitialSolution ("velocity_problem", exactU, 2);
    
    // second step
    exactRho.setTime (t0+dt);
    exactU.setTime (t0+dt);
    exactP.setTime (t0+dt);
    guermondSalgadoMethod.setInitialSolution ("density_problem", exactRho, 1);
    guermondSalgadoMethod.setInitialSolution ("velocity_problem", exactU, 1);
    guermondSalgadoMethod.setInitialSolution ("pressure_update_problem", exactP, 1);

    // define dirichlet boundary conditions
    std::cout << "Define the problem's Dirichlet boundary conditions..." << std::endl;
    dcp::Subdomain dirichletBoundary ((navierstokes::DirichletBoundaryEvaluator ()));
    
    std::cout << "Setting up dirichlet's boundary conditions for 'velocity_problem'" << std::endl;
    guermondSalgadoMethod.addTimeDependentDirichletBC ("velocity_problem", exactU, dirichletBoundary);

    
    // define time dependent external loads
    navierstokes::FirstExternalLoad f1;
    navierstokes::SecondExternalLoad f2;
    guermondSalgadoMethod.problem ("velocity_problem").addTimeDependentCoefficient 
        ("f1", 
         "linear_form", 
         dolfin::reference_to_no_delete_pointer (f1));
    guermondSalgadoMethod.problem ("velocity_problem").addTimeDependentCoefficient 
        ("f2", 
         "linear_form", 
         dolfin::reference_to_no_delete_pointer (f2));
    guermondSalgadoMethod.system ().addLink ("velocity_problem", "rho", "linear_form", "density_problem");
        
    // plots
//    dolfin::plot (mesh, "Mesh");
    
    // solve problem
    std::cout << "Solve the problem..." << std::endl;
    guermondSalgadoMethod.apply ();

    // ------------- //
    // compute error //
    // ------------- //
    
    // computed solution
    auto computedRho = guermondSalgadoMethod.problem ("density_problem").solutions ();
    auto computedU = guermondSalgadoMethod.problem ("velocity_problem").solutions ();
    auto computedP = guermondSalgadoMethod.problem ("pressure_update_problem").solutions ();
    
    // error computers
    error_computers::Form_density_L2_squared_error densityL2SquaredErrorComputer (mesh);
    error_computers::Form_velocity_L2_squared_error velocityL2SquaredErrorComputer (mesh);
    error_computers::Form_velocity_H1_squared_error velocityH1SquaredErrorComputer (mesh);
    error_computers::Form_pressure_L2_squared_error pressureL2SquaredErrorComputer (mesh);
    
    // vectors to contain errors, so that we can get the maximum error later
    std::vector<double> densityL2Errors (computedRho.size ());
    std::vector<double> velocityL2Errors (computedU.size ());
    std::vector<double> velocityH1Errors (computedU.size ());
    std::vector<double> pressureL2Errors (computedP.size ());
    
    std::ofstream ERROR_RHO ("error_rho_test_case_order2_dt=" + std::to_string (dt) + ".txt");
    std::ofstream ERROR_U_H1 ("error_u_H1_test_case_order2_dt=" + std::to_string (dt) + ".txt");
    std::ofstream ERROR_U_L2 ("error_u_L2_test_case_order2_dt=" + std::to_string (dt) + ".txt");
    std::ofstream ERROR_P ("error_p_test_case_order2_dt=" + std::to_string (dt) + ".txt");
    
    std::ofstream MAX_ERROR_RHO ("max_error_rho_test_case_order2.txt", std::ios_base::app);
    std::ofstream MAX_ERROR_U_H1 ("max_error_u_H1_test_case_order2.txt", std::ios_base::app);
    std::ofstream MAX_ERROR_U_L2 ("max_error_u_L2_test_case_order2.txt", std::ios_base::app);
    std::ofstream MAX_ERROR_P ("max_error_p_test_case_order2.txt", std::ios_base::app);
    
    std::cout << "Computing errors..." << std::endl;
    dolfin::Function rescaledComputedP (Q);
    dolfin::Function rescaledPMinusExactP (Q);
    dolfin::Function exactPressure (Q);
    dolfin::Function pressureRescalingFactor (Q);
    
    dolfin::Function compRho (R);
    dolfin::Function compU (V);
    dolfin::Function compP (Q);
    
    for (std::size_t i = 0; i < computedP.size (); ++i)
    {
        const double& rhoTime = computedRho [i+1].first;
        const double& uTime = computedU [i+1].first;
        const double& pTime = computedP [i].first;
            
        compRho = computedRho[i+1].second;
        compU = computedU[i+1].second;
        compP = computedP[i].second;
        
//        dolfin::plot (compRho, "rho, time = " + std::to_string (time));
//        dolfin::plot (compU, "u, time = " + std::to_string (time));
//        dolfin::plot (compP, "p, time = " + std::to_string (time)); 
////        dolfin::plot (computedP[i].second, "p, time = " + std::to_string (time)); 
//        dolfin::interactive ();
//        dolfin::info (computedP[i].second, true);

        exactRho.setTime (rhoTime);
        exactU.setTime (uTime);
        exactP.setTime (pTime);;
        
        densityL2SquaredErrorComputer.exact_rho = exactRho;
        velocityL2SquaredErrorComputer.exact_u = exactU;
        velocityH1SquaredErrorComputer.exact_u = exactU;
        pressureL2SquaredErrorComputer.exact_p = exactP;
    
        densityL2SquaredErrorComputer.rho = computedRho[i].second;
        velocityL2SquaredErrorComputer.u = computedU[i].second;
        velocityH1SquaredErrorComputer.u = computedU[i].second;
    
        dolfin::Array<double> point (2);
        point[0] = 0.0;
        point[1] = 0.0;
        dolfin::Array<double> computedPPointValue (1);
        dolfin::Array<double> exactPPointValue (1);
        computedP[i].second.eval (computedPPointValue, point);
        exactP.eval (exactPPointValue, point);
        double computedPExactPPointDifference = computedPPointValue[0] - exactPPointValue[0];
        pressureRescalingFactor = dolfin::Constant (computedPExactPPointDifference);
        rescaledComputedP = computedP[i].second - pressureRescalingFactor;
//        dolfin::plot (rescaledComputedP, "rescaled computed p, time = " + std::to_string (time));
        exactPressure = exactP;
        rescaledPMinusExactP = rescaledComputedP - exactPressure;
        pressureL2SquaredErrorComputer.p = rescaledComputedP;
        
        densityL2Errors[i] = sqrt (dolfin::assemble (densityL2SquaredErrorComputer));
        velocityL2Errors[i] = sqrt (dolfin::assemble (velocityL2SquaredErrorComputer));
        velocityH1Errors[i] = sqrt (dolfin::assemble (velocityH1SquaredErrorComputer));
        pressureL2Errors[i] = sqrt (dolfin::assemble (pressureL2SquaredErrorComputer));
        
        ERROR_RHO << densityL2Errors[i] << std::endl;
        ERROR_U_H1 << velocityH1Errors[i] << std::endl;
        ERROR_U_L2 << velocityL2Errors[i] << std::endl;
        ERROR_P << pressureL2Errors[i] << std::endl;
    }
    
    std::cout << "done" << std::endl;
    
    double maxDensityL2Error = *(std::max_element (densityL2Errors.begin (), densityL2Errors.end ()));
    double maxVelocityL2Error = *(std::max_element (velocityL2Errors.begin (), velocityL2Errors.end ()));
    double maxVelocityH1Error = *(std::max_element (velocityH1Errors.begin (), velocityH1Errors.end ()));
    double maxPressureL2Error = *(std::max_element (pressureL2Errors.begin (), pressureL2Errors.end ()));
        
    std::cout << "Max density L2 error: " << maxDensityL2Error << std::endl;
    std::cout << "Max velocity L2 error: " << maxVelocityL2Error << std::endl;
    std::cout << "Max velocity H1 error: " << maxVelocityH1Error << std::endl;
    std::cout << "Max pressure L2 error: " << maxPressureL2Error << std::endl;
    
    MAX_ERROR_RHO << dt << " " << maxDensityL2Error << std::endl;
    MAX_ERROR_U_H1 << dt << " " << maxVelocityL2Error << std::endl;
    MAX_ERROR_U_L2 << dt << " " << maxVelocityH1Error << std::endl;
    MAX_ERROR_P << dt << " " << maxPressureL2Error << std::endl;
    // ----------------- //
    // end compute error //
    // ----------------- //
    
//    dolfin::interactive ();
    
    return 0;
}