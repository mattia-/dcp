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
#include "density.h"
#include "velocity.h"
#include "pressure_correction.h"
#include "pressure_update.h"
#include "error_computers.h"
#include "my_div.h"
#include "pressure_correction_2.h"
#include <fstream>

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
    
    class SecondExternalLoadEvaluator
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

double rhoValue = 2.0;
double muValue = 1.0e-1;

namespace exact_solutions
{
    class RhoEvaluator
    {
        public:
            void operator() (dolfin::Array<double>& values, const dolfin::Array<double>& x, const double& t)
            {
                values[0] = rhoValue;
            }
    };
    
    class UEvaluator
    {
        public:
            void operator() (dolfin::Array<double>& values, const dolfin::Array<double>& x, const double& t)
            {
                values[0] = x[1] * (1 - x[1]);
                values[1] = 0;
            }
    };
    
    class PEvaluator
    {
        public:
            void operator() (dolfin::Array<double>& values, const dolfin::Array<double>& x, const double& t)
            {
                values[0] = -2 * muValue / rhoValue *x[0] + 2 * muValue / rhoValue * 5;
            }
    };
}

int main (int argc, char* argv[])
{
//    dolfin::set_log_level (dolfin::DBG);
    // create mesh and finite element space
    std::cout << "Create mesh and finite element space..." << std::endl;
    dolfin::RectangleMesh mesh (0.0, 0.0, 5.0, 1.0, 40, 40);
    
    density::FunctionSpace R (mesh);
    velocity::FunctionSpace V (mesh);
    pressure_correction::FunctionSpace Q (mesh);
    
    // define constant
    std::cout << "Define the coefficients..." << std::endl;
    dolfin::Constant mu (muValue);
    dolfin::Constant chi (1.0);
    double t0 = 0.0;
    double dt = 0.1;
    dolfin::Constant dtConst (dt);
    double T = 1.0;
    
    // exact solution
    dcp::TimeDependentExpression exactRho ((exact_solutions::RhoEvaluator ()));
    exactRho.setTime (t0);
    dcp::TimeDependentExpression exactU (2, exact_solutions::UEvaluator ());
    exactU.setTime (t0);
    dcp::TimeDependentExpression exactP ((exact_solutions::PEvaluator ()));
    exactP.setTime (t0);

    // define the problems
    std::cout << "Define the problem..." << std::endl;
    dcp::GuermondSalgadoMethod <density::BilinearForm,
                                density::LinearForm,
                                velocity::BilinearForm,
                                velocity::LinearForm,
                                pressure_correction::BilinearForm,
                                pressure_correction::LinearForm,
                                pressure_update::BilinearForm,
                                pressure_update::LinearForm>
        guermondSalgadoMethod (R, V, Q, t0, dt, T, mu, chi);
    
    // set initial solutions
    guermondSalgadoMethod.setInitialSolution ("density_problem", exactRho);
    guermondSalgadoMethod.setInitialSolution ("velocity_problem", exactU);
    guermondSalgadoMethod.setInitialSolution ("pressure_update_problem", exactP);
    
    
    dolfin::plot (guermondSalgadoMethod.problem ("pressure_update_problem").solution (), "p iniziale");
    dolfin::plot (guermondSalgadoMethod.problem ("velocity_problem").solution (), "u iniziale");
//    dolfin::interactive ();
    
    // ******************************************************** //
    guermondSalgadoMethod.system ().removeLink ("velocity_problem", "rho", "bilinear_form");
    guermondSalgadoMethod.system ().removeLinkToPreviousSolution ("velocity_problem", "rho_old", "bilinear_form"); 
    guermondSalgadoMethod.system ().removeLinkToPreviousSolution ("velocity_problem", "rho_old", "linear_form"); 
    guermondSalgadoMethod.problem ("velocity_problem").setCoefficient ("bilinear_form", dolfin::reference_to_no_delete_pointer (exactRho),  "rho");
    guermondSalgadoMethod.problem ("velocity_problem").setCoefficient ("bilinear_form", dolfin::reference_to_no_delete_pointer (exactRho),  "rho_old");
    guermondSalgadoMethod.problem ("velocity_problem").setCoefficient ("linear_form", dolfin::reference_to_no_delete_pointer (exactRho),  "rho_old");
    // ******************************************************** //
    
    // define dirichlet boundary conditions
    std::cout << "Define the problem's Dirichlet boundary conditions..." << std::endl;
    dcp::Subdomain noSlipBoundary ((navierstokes::NoSlipBoundaryEvaluator ()));
    dcp::Subdomain inflowBoundary ((navierstokes::InflowBoundaryEvaluator ()));
    dcp::Subdomain outflowBoundary ((navierstokes::OutflowBoundaryEvaluator ()));
    dolfin::Constant zero (0.0);
    
    std::cout << "Setting up dirichlet's boundary conditions for 'velocity_problem'" << std::endl;
    guermondSalgadoMethod.addDirichletBC ("density_problem", exactRho, inflowBoundary);
    guermondSalgadoMethod.addDirichletBC ("velocity_problem", exactU, noSlipBoundary);
    guermondSalgadoMethod.addDirichletBC ("velocity_problem", exactU, inflowBoundary);
    guermondSalgadoMethod.addDirichletBC ("pressure_correction_problem", zero, outflowBoundary);

    dcp::TimeDependentExpression f1 (2, (navierstokes::ExternalLoadEvaluator ()));
    dcp::TimeDependentExpression f2 (2, (navierstokes::SecondExternalLoadEvaluator ()));
    guermondSalgadoMethod.problem ("velocity_problem").addTimeDependentCoefficient 
        ("f1", 
         "linear_form", 
         dolfin::reference_to_no_delete_pointer (f1));
    guermondSalgadoMethod.problem ("velocity_problem").addTimeDependentCoefficient 
        ("f2", 
         "linear_form", 
         dolfin::reference_to_no_delete_pointer (f2));
    guermondSalgadoMethod.system ().addLink ("velocity_problem", "rho", "linear_form", "density_problem");
        
    
    // define time dependent external loads
    for (auto& name : guermondSalgadoMethod.system ().problemsNames ())
    {
//        guermondSalgadoMethod.problem (name).parameters ["plot_interval"] = 1;
//        guermondSalgadoMethod.problem (name).parameters ["plot_title"] = name;
    }
//    guermondSalgadoMethod.problem ("pressure_update_problem").parameters ["pause"] = true;
//        guermondSalgadoMethod.problem ("velocity_problem").parameters ["plot_interval"] = 1;
//        guermondSalgadoMethod.problem ("velocity_problem").parameters ["plot_title"] = "velocity_problem";
//        guermondSalgadoMethod.problem ("pressure_update_problem").parameters ["plot_interval"] = 1;
//        guermondSalgadoMethod.problem ("pressure_update_problem").parameters ["plot_title"] = "pressure_update_problem";
//        guermondSalgadoMethod.problem ("pressure_correction_problem").parameters ["plot_interval"] = 1;
//        guermondSalgadoMethod.problem ("pressure_correction_problem").parameters ["plot_title"] = "pressure_correction_problem";
    
    // plots
    dolfin::plot (mesh, "Mesh");
    
    // solve problem
    std::cout << "Solve the problem..." << std::endl;
    guermondSalgadoMethod.apply ();

    // ------------- //
    // compute error //
    // ------------- //
    
    // computed solution
    auto computedRho = guermondSalgadoMethod.problem ("density_problem").solutionsWithTimes ();
    auto computedU = guermondSalgadoMethod.problem ("velocity_problem").solutionsWithTimes ();
    auto computedPhi = guermondSalgadoMethod.problem ("pressure_correction_problem").solutionsWithTimes ();
    auto computedP = guermondSalgadoMethod.problem ("pressure_update_problem").solutionsWithTimes ();
    
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
    
    std::ofstream ERROR_U_H1 ("errore_u_H1.txt");
    std::ofstream ERROR_U_L2 ("errore_u_L2.txt");
    std::ofstream ERROR_P ("errore_p.txt");
    
    std::cout << "Computing errors..." << std::endl;
    dolfin::Function rescaledComputedP (Q);
    dolfin::Function rescaledPDifference (Q);
    dolfin::Function exactPressure (Q);
    dolfin::Function pressurePointDifference (Q);
    
    my_div::BilinearForm a (Q, Q);
    my_div::LinearForm L (Q);
    dolfin::Function sol (Q);
    
    pressure_correction_2::BilinearForm aPressureCorrection2 (Q, Q);
    pressure_correction_2::LinearForm LPressureCorrection2 (Q);
    dolfin::Function pressureCorrection2 (Q);
    
    dolfin::Function compRho (R);
    dolfin::Function compU (V);
    dolfin::Function compPhi (Q);
    dolfin::Function compP (Q);
    
    dolfin::Function phiDiff (Q);
    dolfin::Function rescalingForDivU (Q);
    dolfin::Function rescaledDivU (Q);
    
    for (std::size_t i = 0; i < computedRho.size (); ++i)
    {
        const double& time = computedRho [i].first;
        std::cout << "****************** t = " << time << std::endl; 
            
        compRho = *(computedRho[i].second);
        compU = *(computedU[i].second);
        compPhi = *(computedPhi[i].second);
        compP = *(computedP[i].second);
        dolfin::plot (compRho, "rho, time = " + std::to_string (time));
        dolfin::plot (compU, "u, time = " + std::to_string (time));
        dolfin::plot (compPhi, "phi, time = " + std::to_string (time)); 
        dolfin::plot (compP, "p, time = " + std::to_string (time)); 


        exactRho.setTime (time);
        exactU.setTime (time);
        exactP.setTime (time);
        
        densityL2SquaredErrorComputer.exact_rho = exactRho;
        velocityL2SquaredErrorComputer.exact_u = exactU;
        velocityH1SquaredErrorComputer.exact_u = exactU;
        pressureL2SquaredErrorComputer.exact_p = exactP;
    
        densityL2SquaredErrorComputer.rho = *(computedRho[i].second);
        velocityL2SquaredErrorComputer.u = *(computedU[i].second);
        velocityH1SquaredErrorComputer.u = *(computedU[i].second);
        pressureL2SquaredErrorComputer.p = *(computedP[i].second);
        
//        dolfin::Array<double> point (2);
//        point[0] = 0.0;
//        point[1] = 0.0;
//        dolfin::Array<double> computedPPointValue (1);
//        dolfin::Array<double> exactPPointValue (1);
//        computedP[i].second->eval (computedPPointValue, point);
//        exactP.eval (exactPPointValue, point);
//        std::cout << "exactPPointValue = " << exactPPointValue[0] << std::endl;
//        std::cout << "computedPPointValue = " << computedPPointValue[0] << std::endl;
//        double difference = computedPPointValue[0] - exactPPointValue[0];
//        std::cout << "pressurePointDifference = " << difference << std::endl;
//        dolfin::Constant computedExactPPointDifference (difference);
//        pressurePointDifference = computedExactPPointDifference;
//        dolfin::plot (pressurePointDifference, "pressurePointDifference");
//        std::cout << "rescaledComputedP = " << computedPPointValue[0] - difference << std::endl;
//        rescaledComputedP = *(computedP[i].second) - pressurePointDifference;
//        dolfin::plot (rescaledComputedP, "rescaledComputedP");
//        exactPressure = exactP;
//        rescaledPDifference = rescaledComputedP - exactPressure;
//        dolfin::plot (rescaledPDifference, "rescaledPDifference");
//        pressureL2SquaredErrorComputer.p = rescaledComputedP;
        
    
        densityL2Errors[i] = sqrt (dolfin::assemble (densityL2SquaredErrorComputer));
        velocityL2Errors[i] = sqrt (dolfin::assemble (velocityL2SquaredErrorComputer));
        velocityH1Errors[i] = sqrt (dolfin::assemble (velocityH1SquaredErrorComputer));
        pressureL2Errors[i] = sqrt (dolfin::assemble (pressureL2SquaredErrorComputer));
        ERROR_U_H1 << velocityH1Errors[i] << std::endl;
        ERROR_U_L2 << velocityL2Errors[i] << std::endl;
        ERROR_P << pressureL2Errors[i] << std::endl;
        
        L.set_coefficient ("vel", computedU[i].second);
        solve (a == L, sol);
        dolfin::plot (sol, "div(u), time = " + std::to_string(time));
        
//        double rescalingFactor;
//        std::cout << "rescaling factor for div u =  ";
//        std::cin >> rescalingFactor;
//        dolfin::Constant rescalingFactorAsDolfinConstant (rescalingFactor);
//        rescalingForDivU = rescalingFactorAsDolfinConstant;
//        dolfin::plot (rescalingForDivU, "rescalingForDivU");
//        
//        rescaledDivU = sol - rescalingForDivU;
        rescaledDivU = sol;
//        dolfin::plot (rescaledDivU, "div(u) rescaled, time = " + std::to_string(time));
        LPressureCorrection2.set_coefficient ("divu", dolfin::reference_to_no_delete_pointer (rescaledDivU));
        LPressureCorrection2.set_coefficient ("chi", dolfin::reference_to_no_delete_pointer (chi));
        LPressureCorrection2.set_coefficient ("dt", dolfin::reference_to_no_delete_pointer (dtConst));
        dolfin::DirichletBC bc (Q, zero, outflowBoundary);
        solve (aPressureCorrection2 == LPressureCorrection2, pressureCorrection2, bc);
//        dolfin::plot (pressureCorrection2, "pressure correction");
                
        phiDiff = pressureCorrection2 - compPhi;
//        dolfin::plot (phiDiff, "phi diff");
               
        dolfin::interactive ();
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
    // ----------------- //
    // end compute error //
    // ----------------- //
    
    dolfin::interactive ();
    
    return 0;
}
