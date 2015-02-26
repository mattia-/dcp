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
#include "chorin_temam_prediction.h"
#include "chorin_temam_correction.h"
#include "chorin_temam_velocity_projection.h"
#include "chorin_temam_pressure_update.h"
#include "error_computers.h"
#include <fstream>

double rhoValue = 2.0;
double muValue = 1.0e-1;

namespace navierstokes
{
    class ExternalLoadEvaluator
    {
        public:
            void operator() (dolfin::Array<double>& values, const dolfin::Array<double>& x, const double& t)
            {
                values[0] = 0;
                values[1] = 0;
            }
    }; 
    
    class NoSlipBoundaryEvaluator
    { 
        public:
            bool operator() (const dolfin::Array<double>& x, bool on_boundary)
            {
                return (dolfin::near (x[1], 0) && on_boundary) || (dolfin::near (x[1], 1) && on_boundary);
            }
    }; 
    
    class InflowBoundaryEvaluator
    { 
        public:
            bool operator() (const dolfin::Array<double>& x, bool on_boundary)
            {
                return (dolfin::near (x[0], 0) && on_boundary);
            }
    }; 
    
    class OutflowBoundaryEvaluator
    { 
        public:
            bool operator() (const dolfin::Array<double>& x, bool on_boundary)
            {
                return (dolfin::near (x[0], 5) && on_boundary);
            }
    }; 
}

namespace exact_solutions
{
    class UEvaluator
    {
        public:
            void operator() (dolfin::Array<double>& values, const dolfin::Array<double>& x)
            {
                values[0] = x[1] * (1 - x[1]);
                values[1] = 0;
            }
    };
    
    class PEvaluator
    {
        public:
            void operator() (dolfin::Array<double>& values, const dolfin::Array<double>& x)
            {
                values[0] = -2 * muValue / rhoValue * x[0] + 2 * muValue / rhoValue * 5;
            }
    };
}

int main (int argc, char* argv[])
{
//    dolfin::set_log_level (dolfin::DBG);
    // create mesh and finite element space
    std::cout << "Create mesh and finite element space..." << std::endl;
    dolfin::RectangleMesh mesh (0.0, 0.0, 5.0, 1.0, 40, 40);
    
    chorin_temam_prediction::FunctionSpace V (mesh);
    chorin_temam_correction::FunctionSpace Q (mesh);
    
    // define constant
    std::cout << "Define the coefficients..." << std::endl;
    dolfin::Constant nu (muValue / rhoValue);
    double t0 = 0.0;
    double dt = 0.1;
    double T = 1.0;
    
    // exact solution
    dcp::Expression exactU (2, exact_solutions::UEvaluator ());
    dcp::Expression exactP ((exact_solutions::PEvaluator ()));

    // define the problems
    std::cout << "Define the problem..." << std::endl;
    dcp::IncrementalChorinTemamMethod <chorin_temam_prediction::BilinearForm,
                                       chorin_temam_prediction::LinearForm,
                                       chorin_temam_correction::BilinearForm,
                                       chorin_temam_correction::LinearForm,
                                       chorin_temam_velocity_projection::BilinearForm,
                                       chorin_temam_velocity_projection::LinearForm,
                                       chorin_temam_pressure_update::BilinearForm,
                                       chorin_temam_pressure_update::LinearForm>
        incrementalChorinTemamMethod (V, Q, t0, dt, T, nu);
    
    // set initial solutions
    incrementalChorinTemamMethod.setInitialSolution ("prediction_problem", exactU);
    incrementalChorinTemamMethod.setInitialSolution ("correction_problem", exactP);

    dolfin::Constant zero (0.0);
    
    // define dirichlet boundary conditions
    std::cout << "Define the problem's Dirichlet boundary conditions..." << std::endl;
    dcp::Subdomain noSlipBoundary ((navierstokes::NoSlipBoundaryEvaluator ()));
    dcp::Subdomain inflowBoundary ((navierstokes::InflowBoundaryEvaluator ()));
    dcp::Subdomain outflowBoundary ((navierstokes::OutflowBoundaryEvaluator ()));
    
    std::cout << "Setting up dirichlet's boundary conditions" << std::endl;
    incrementalChorinTemamMethod.addDirichletBC ("prediction_problem", exactU, noSlipBoundary);
    incrementalChorinTemamMethod.addDirichletBC ("prediction_problem", exactU, inflowBoundary);
    incrementalChorinTemamMethod.addDirichletBC ("correction_problem", zero, outflowBoundary);

    dcp::TimeDependentExpression f (2, (navierstokes::ExternalLoadEvaluator ()));
    incrementalChorinTemamMethod.problem ("prediction_problem").addTimeDependentCoefficient 
        ("f", 
         "linear_form", 
         dolfin::reference_to_no_delete_pointer (f));
        
    // plots
    dolfin::plot (mesh, "Mesh");
    
    // solve problem
    std::cout << "Solve the problem..." << std::endl;
    incrementalChorinTemamMethod.apply ();

    // ------------- //
    // compute error //
    // ------------- //
    
    // computed solution
    auto computedU = incrementalChorinTemamMethod.problem ("velocity_projection_problem").solutions ();
    auto computedP = incrementalChorinTemamMethod.problem ("pressure_update_problem").solutions ();
    
    // error computers
    error_computers::Form_velocity_L2_squared_error velocityL2SquaredErrorComputer (mesh);
    error_computers::Form_velocity_H1_squared_error velocityH1SquaredErrorComputer (mesh);
    error_computers::Form_pressure_L2_squared_error pressureL2SquaredErrorComputer (mesh);
    
    // vectors to contain errors, so that we can get the maximum error later
    std::vector<double> velocityL2Errors (computedU.size ());
    std::vector<double> velocityH1Errors (computedU.size ());
    std::vector<double> pressureL2Errors (computedP.size ());
    
    std::ofstream ERROR_U_H1 ("errore_u_H1_poiseuille_chorin_temam.txt");
    std::ofstream ERROR_U_L2 ("errore_u_L2_poiseuille_chorin_temam.txt");
    std::ofstream ERROR_P ("errore_p_poiseuille_chorin_temam.txt");
    
    dolfin::Function compU (V);
    dolfin::Function compP (Q);
    
    std::cout << "Computing errors..." << std::endl;
    
    for (std::size_t i = 0; i < computedU.size (); ++i)
    {
        const double& time = computedU [i].first;
            
        compU = computedU[i].second;
        compP = computedP[i].second;
        dolfin::plot (compU, "u, time = " + std::to_string (time));
        dolfin::plot (compP, "p, time = " + std::to_string (time)); 
        
        velocityL2SquaredErrorComputer.exact_u = exactU;
        velocityH1SquaredErrorComputer.exact_u = exactU;
        pressureL2SquaredErrorComputer.exact_p = exactP;
    
        velocityL2SquaredErrorComputer.u = computedU[i].second;
        velocityH1SquaredErrorComputer.u = computedU[i].second;
        pressureL2SquaredErrorComputer.p = computedP[i].second;
        
        velocityL2Errors[i] = sqrt (dolfin::assemble (velocityL2SquaredErrorComputer));
        velocityH1Errors[i] = sqrt (dolfin::assemble (velocityH1SquaredErrorComputer));
        pressureL2Errors[i] = sqrt (dolfin::assemble (pressureL2SquaredErrorComputer));
        
        ERROR_U_H1 << velocityH1Errors[i] << std::endl;
        ERROR_U_L2 << velocityL2Errors[i] << std::endl;
        ERROR_P << pressureL2Errors[i] << std::endl;
        
//        dolfin::interactive ();
    }
    
    std::cout << "done" << std::endl;
    
    double maxVelocityL2Error = *(std::max_element (velocityL2Errors.begin (), velocityL2Errors.end ()));
    double maxVelocityH1Error = *(std::max_element (velocityH1Errors.begin (), velocityH1Errors.end ()));
    double maxPressureL2Error = *(std::max_element (pressureL2Errors.begin (), pressureL2Errors.end ()));
        
    std::cout << "Max velocity L2 error: " << maxVelocityL2Error << std::endl;
    std::cout << "Max velocity H1 error: " << maxVelocityH1Error << std::endl;
    std::cout << "Max pressure L2 error: " << maxPressureL2Error << std::endl;
    // ----------------- //
    // end compute error //
    // ----------------- //
    
    dolfin::interactive ();
    
    return 0;
}
