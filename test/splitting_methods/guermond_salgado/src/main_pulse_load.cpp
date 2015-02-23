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
#include <fstream>

double rhoValue = 1.0;
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
    
    class SecondExternalLoadEvaluator
    {
        public:
            void operator() (dolfin::Array<double>& values, const dolfin::Array<double>& x, const double& t)
            {
                values[0] = -cos (x[0]) * sin (x[1]) * (2 * cos (2*t) + 2 * muValue / rhoValue * sin (2*t)); 
                values[1] =  sin (x[0]) * cos (x[1]) * (2 * cos (2*t) + 2 * muValue / rhoValue * sin (2*t)); 
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
                values[0] = rhoValue;
            }
    };
    
    class UEvaluator
    {
        public:
            void operator() (dolfin::Array<double>& values, const dolfin::Array<double>& x, const double& t)
            {
                values[0] = - cos (x[0]) * sin (x[1]) * sin (2*t);
                values[1] =   sin (x[0]) * cos (x[1]) * sin (2*t);
            }
    };
    
    class PEvaluator
    {
        public:
            void operator() (dolfin::Array<double>& values, const dolfin::Array<double>& x, const double& t)
            {
                values[0] = -1.0/4.0 * (cos (2*x[0]) + cos (2*x[1])) * sin (2*t) * sin (2*t);
            }
    };
}

int main (int argc, char* argv[])
{
//    dolfin::set_log_level (dolfin::DBG);
    // create mesh and finite element space
    std::cout << "Create mesh and finite element space..." << std::endl;
//    mshr::Circle circle (dolfin::Point (0.0, 0.0), 1);
//    dolfin::Mesh mesh;
//    mshr::generate (mesh, circle, 50);
    dolfin::UnitSquareMesh mesh (40, 40);
    
    density::FunctionSpace R (mesh);
    velocity::FunctionSpace V (mesh);
    pressure_correction::FunctionSpace Q (mesh);
    
    // define constant
    std::cout << "Define the coefficients..." << std::endl;
    dolfin::Constant mu (muValue);
    dolfin::Constant chi (1.0);
    double t0 = 0.0;
    double dt = 0.1;
    double T = 1;
    
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
    dcp::Subdomain noSlipBoundary ((navierstokes::DirichletBoundaryEvaluator ()));
    
    std::cout << "Setting up dirichlet's boundary conditions for 'velocity_problem'" << std::endl;
    guermondSalgadoMethod.addTimeDependentDirichletBC ("velocity_problem", exactU, noSlipBoundary);

    
    // define time dependent external loads
    dcp::TimeDependentExpression f1 (2, (navierstokes::ExternalLoadEvaluator ()));
    dcp::TimeDependentExpression f2 (2, (navierstokes::SecondExternalLoadEvaluator ()));
    dolfin::Constant zero (0.0);
    guermondSalgadoMethod.problem ("velocity_problem").addTimeDependentCoefficient 
        ("f1", 
         "linear_form", 
         dolfin::reference_to_no_delete_pointer (f1));
    guermondSalgadoMethod.problem ("velocity_problem").addTimeDependentCoefficient 
        ("f2", 
         "linear_form", 
         dolfin::reference_to_no_delete_pointer (f2));
    guermondSalgadoMethod.problem ("velocity_problem").setCoefficient 
        ("linear_form", 
         dolfin::reference_to_no_delete_pointer (zero),
         "rho");
        
    
    for (auto& name : guermondSalgadoMethod.system ().problemsNames ())
    {
//        guermondSalgadoMethod.problem (name).parameters ["plot_interval"] = 1;
//        guermondSalgadoMethod.problem (name).parameters ["plot_title"] = name;
    }
//    guermondSalgadoMethod.problem ("pressure_update_problem").parameters ["pause"] = true;
    
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
    std::size_t i;
    for (i = 0; i < computedRho.size (); ++i)
    {
        const double& time = computedRho [i].first;
        std::cout << "****************** t = " << time << std::endl; 
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
        dolfin::Array<double> point (2);
        point[0] = 0.0;
        point[1] = 0.0;
        dolfin::Array<double> computedPPointValue (1);
        dolfin::Array<double> exactPPointValue (1);
        computedP[i].second->eval (computedPPointValue, point);
        exactP.eval (exactPPointValue, point);
        std::cout << "exactPPointValue = " << exactPPointValue[0] << std::endl;
        std::cout << "computedPPointValue = " << computedPPointValue[0] << std::endl;
        double difference = computedPPointValue[0] - exactPPointValue[0];
        std::cout << "pressurePointDifference = " << difference << std::endl;
        dolfin::Constant computedExactPPointDifference (difference);
        pressurePointDifference = computedExactPPointDifference;
        dolfin::plot (pressurePointDifference, "pressurePointDifference");
        std::cout << "rescaledComputedP = " << computedPPointValue[0] - difference << std::endl;
        rescaledComputedP = *(computedP[i].second) - pressurePointDifference;
        dolfin::plot (rescaledComputedP, "rescaledComputedP");
        exactPressure = exactP;
        rescaledPDifference = rescaledComputedP - exactPressure;
        dolfin::plot (rescaledPDifference, "rescaledPDifference");
        pressureL2SquaredErrorComputer.p = rescaledComputedP;
        
//        dolfin::interactive ();
    
        densityL2Errors[i] = sqrt (dolfin::assemble (densityL2SquaredErrorComputer));
        velocityL2Errors[i] = sqrt (dolfin::assemble (velocityL2SquaredErrorComputer));
        velocityH1Errors[i] = sqrt (dolfin::assemble (velocityH1SquaredErrorComputer));
        pressureL2Errors[i] = sqrt (dolfin::assemble (pressureL2SquaredErrorComputer));
        ERROR_U_H1 << velocityH1Errors[i] << std::endl;
        ERROR_U_L2 << velocityL2Errors[i] << std::endl;
        ERROR_P << pressureL2Errors[i] << std::endl;
    }
    
    dolfin::Function uDifference (V);
    dolfin::Function uEsatta (V);
    uEsatta = exactU;
    uDifference = *(computedU[i-1].second) - uEsatta;
    dolfin::plot (uDifference, "differenza u");
    
    dolfin::Function pDifference (Q);
    dolfin::Function pEsatta (Q);
    pEsatta = exactP;
    pDifference = *(computedP[i-1].second) - pEsatta;
    dolfin::plot (pDifference, "differenza p");
    
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
