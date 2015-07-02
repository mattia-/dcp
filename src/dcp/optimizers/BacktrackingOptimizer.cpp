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

#include <dcp/optimizers/BacktrackingOptimizer.h>
#include <dolfin/parameter/Parameters.h>
#include <dolfin/log/dolfin_log.h>
#include <dolfin/function/Function.h>
#include <dolfin/common/NoDeleter.h>
#include <dolfin/fem/assemble.h>
#include <dolfin/fem/Form.h>
#include <string>
#include <functional>
#include <cmath>
#include <iomanip>
#include <dcp/utils/dotproductforms.h>

namespace dcp
{
    BacktrackingOptimizer::BacktrackingOptimizer ():
        GenericDescentMethod ()
    {
        dolfin::begin (dolfin::DBG, "Creating BacktrackingOptimizer object...");
        
        dolfin::log (dolfin::DBG, "Setting up parameters...");
        parameters.add ("descent_method", "backtracking_gradient_method");
        parameters.add ("gradient_norm_tolerance", 1e-6);
        parameters.add ("relative_increment_tolerance", 1e-6);
        parameters.add ("convergence_criterion", "both");
        parameters.add ("c_1", 1e-3);
        parameters.add ("alpha_0", 1.0);
        parameters.add ("rho", 0.5);
        parameters.add ("max_minimization_iterations", 100);
        parameters.add ("max_backtracking_iterations", 20);
        parameters.add ("output_file_name", "");
        
        dolfin::log (dolfin::DBG, "BacktrackingOptimizer object created");
        
        dolfin::end ();
    }
    


    void BacktrackingOptimizer::apply (dcp::GenericEquationSystem& problem,
                                       const dcp::GenericObjectiveFunctional& objectiveFunctional, 
                                       dolfin::Function& initialGuess,
                                       const dcp::GenericDescentMethod::Updater& updater)
    {
        // ------------------------------------------------------------------------------------------------------- //
        // minimization loop settings
        // ------------------------------------------------------------------------------------------------------- //
        // get parameters values
        double c_1                        = this->parameters ["c_1"];
        double alpha_0                    = this->parameters ["alpha_0"];
        double rho                        = this->parameters ["rho"];
        double gradientNormTolerance      = this->parameters ["gradient_norm_tolerance"];
        double relativeIncrementTolerance = this->parameters ["relative_increment_tolerance"];
        std::string convergenceCriterion  = this->parameters ["convergence_criterion"];
        int maxMinimizationIterations     = this->parameters ["max_minimization_iterations"];
        int maxBacktrackingIterations     = this->parameters ["max_backtracking_iterations"];

        // define loop variables
        double alpha;
        double gradientNorm;
        double gradientDotSearchDirection;
        double incrementNorm;
        double previousControlVariableNorm;
        double relativeIncrement;
        double previousFunctionalValue;
        double currentFunctionalValue;
        dolfin::Function& controlVariable = initialGuess;
        dolfin::Function functionalGradient (controlVariable.function_space ());
        dolfin::Function searchDirection (controlVariable.function_space ());
        dolfin::Function controlVariableIncrement (controlVariable.function_space ());
        
        // define output file and print header if necessary
        std::ofstream OUTFILE;
        bool hasOutputFile = openOutputFile (OUTFILE);
        
        
        // all linear problems in the EquationSystem "problem" should be reassembled every time. So we
        // set the parameter "force_reassemble_system" to true for every one of them
        dolfin::begin (dolfin::DBG, "Scanning composite differential problem...");
        for (std::size_t i = 0; i < problem.size (); ++i)
        {
            if (static_cast<std::string> (problem[i].parameters["problem_type"]) == "linear")
            {
                dolfin::log (dolfin::DBG, "Set parameter \"force_reassemble_system\" to true for problem number %d", i);
                problem[i].parameters["force_reassemble_system"] = true;
            }
        }
        dolfin::end ();
        
        
        // function to check convergence
        std::function <bool ()> isConverged;
        
        // set value of isConverged () according to convergence criterion
        if (convergenceCriterion == "gradient")
        {
            isConverged = [&] () {return gradientNorm < gradientNormTolerance;};
        }
        else if (convergenceCriterion == "increment")
        {
            isConverged = [&] () {return relativeIncrement < relativeIncrementTolerance;};
        }
        else if (convergenceCriterion == "both")
        {
            isConverged = [&] () {return (relativeIncrement < relativeIncrementTolerance) || (gradientNorm < gradientNormTolerance);};
        }
        else
        {
            dolfin::dolfin_error ("dcp: BacktrackingOptimizer.cpp",
                                  "apply",
                                  "Unknown convergence criterion \"%s\"", 
                                  convergenceCriterion.c_str ());
        }
        // ------------------------------------------------------------------------------------------------------- //
        // end of minimization loop setting
        // ------------------------------------------------------------------------------------------------------- //

        
        // update and solve the problem for the first time
        dolfin::begin ("Minimization loop initialization...");
        
        dolfin::begin (dolfin::PROGRESS, "Updating differential problem using control initial guess...");
        updater (problem, controlVariable);
        dolfin::end ();

        // solve problems
        problem.solve ();
        
        
        // initialize loop variables
        dolfin::begin (dolfin::PROGRESS, "Computing norm of functional gradient...");
        gradientNorm = dotProduct_.norm (objectiveFunctional.gradient (), 
                                         *(controlVariable.function_space ()->mesh ()));
        dolfin::end ();
        
        relativeIncrement = relativeIncrementTolerance + 1; // just an initialization to be sure that the first iteration 
                                                            // of the minimization loop is performed
        dolfin::begin (dolfin::PROGRESS, "Evaluating functional...");
        currentFunctionalValue = objectiveFunctional.evaluateFunctional ();
        dolfin::end ();

        dolfin::end (); // "Minimization loop initialization"

        dolfin::log (dolfin::DBG, "=================\nMINIMIZATION LOOP\n=================");
        dolfin::log (dolfin::DBG, "Parameters value are:");
        dolfin::log (dolfin::DBG, "c_1 = %f", c_1);
        dolfin::log (dolfin::DBG, "alpha_0 = %f", alpha_0);
        dolfin::log (dolfin::DBG, "rho = %f", rho);
        dolfin::log (dolfin::DBG, "gradient norm tolerance = %f", gradientNormTolerance);
        dolfin::log (dolfin::DBG, "relative increment tolerance = %f", relativeIncrementTolerance);
        dolfin::log (dolfin::DBG, "convergence criterion = %s", convergenceCriterion.c_str ());
        dolfin::log (dolfin::DBG, "maximum minimization iterations = %d", maxMinimizationIterations);
        dolfin::log (dolfin::DBG, "maximum backtracking iterations = %d", maxBacktrackingIterations);
        

        dolfin::log (dolfin::INFO, "");
        dolfin::log (dolfin::INFO, 
                     "============================\nLOOP VARIABLES INITIAL VALUE\n============================");
        dolfin::log (dolfin::INFO, "gradient norm = %f", gradientNorm);
        dolfin::log (dolfin::INFO, "functional value = %f\n", currentFunctionalValue);
            
        dolfin::log (dolfin::INFO, "***************************"); 
        dolfin::log (dolfin::INFO, "**** MINIMIZATION LOOP ****");
        dolfin::begin ("***************************"); 
        int minimizationIteration = 0;
        
        // print results to file
        if (hasOutputFile)
        {
            print (OUTFILE, minimizationIteration, currentFunctionalValue, 0, 0, gradientNorm, relativeIncrement);
        }
        
        while (isConverged () == false && minimizationIteration < maxMinimizationIterations)
        {
            minimizationIteration++;
            
            dolfin::log (dolfin::INFO, "==========================");
            dolfin::log (dolfin::INFO, "Minimization iteration %d", minimizationIteration);
            dolfin::begin ("==========================");
            
            // iteration-specific variable initialization
            alpha = alpha_0;
            int backtrackingIteration = 0;
            previousFunctionalValue = currentFunctionalValue;
            functionalGradient = objectiveFunctional.gradient ();
            dolfin::begin (dolfin::PROGRESS, "Computing search direction...");
            searchDirectionComputer_ (searchDirection, functionalGradient);
            dolfin::end ();
            
            // compute dot product between gradient and search direction
            dolfin::begin (dolfin::PROGRESS, "Computing dot product between gradient and search direction...");
            gradientDotSearchDirection = dotProduct_.compute (objectiveFunctional.gradient (), 
                                                              searchDirection,
                                                              *(controlVariable.function_space ()->mesh ()));
            dolfin::end ();
            
            // solution of problem with alpha_0
            dolfin::log (dolfin::INFO, "Alpha = %f", alpha);

            // update control variable value
            dolfin::begin (dolfin::PROGRESS, "Updating control variable...");
            dolfin::Function previousControlVariable (controlVariable);
            controlVariable = previousControlVariable + (searchDirection * alpha);
            dolfin::end ();

            // update problem 
            dolfin::begin (dolfin::PROGRESS, "Updating differential problem...");
            updater (problem, controlVariable);
            dolfin::end ();

            // solve problem
            problem.solve ();

            // update value of the functional
            dolfin::begin (dolfin::PROGRESS, "Evaluating functional...");
            currentFunctionalValue = objectiveFunctional.evaluateFunctional ();
            dolfin::end ();
            
            dolfin::log (dolfin::INFO, "Functional value = %f\n", currentFunctionalValue);
            
            
            // backtracking loop
            backtrackingLoop (previousFunctionalValue, 
                              currentFunctionalValue, 
                              gradientDotSearchDirection, 
                              alpha,
                              backtrackingIteration,
                              controlVariable,
                              previousControlVariable,
                              searchDirection,
                              problem,
                              objectiveFunctional,
                              updater);
             
            // update gradient norm if necessary for convergence check
            if (convergenceCriterion == "gradient" || convergenceCriterion == "both")
            {
                dolfin::begin (dolfin::PROGRESS, "Computing norm of functional gradient...");
                gradientNorm = dotProduct_.norm (objectiveFunctional.gradient (), 
                                                 *(controlVariable.function_space ()->mesh ()));
                dolfin::end ();
                
                dolfin::log (dolfin::INFO, "Gradient norm = %f", gradientNorm);
            }
            

            // update relative increment if necessary for convergence check
            if (convergenceCriterion == "increment" || convergenceCriterion == "both")
            {
                dolfin::begin (dolfin::PROGRESS, "Computing relative increment...");
                previousControlVariableNorm = dotProduct_.norm (previousControlVariable);
                dolfin::end ();

                controlVariableIncrement = controlVariable - previousControlVariable;
                incrementNorm = dotProduct_.norm (controlVariableIncrement);

                // compute relative increment. We add DOLFIN_EPS at the denominator in case previousControlVariableNorm 
                // is zero (like at the first iteration)
                relativeIncrement = incrementNorm / (previousControlVariableNorm + DOLFIN_EPS); 
                dolfin::log (dolfin::INFO, "Relative increment = %f", relativeIncrement);
            }
                
            dolfin::log (dolfin::INFO, "Functional value = %f\n", currentFunctionalValue);
            
            if (hasOutputFile)
            {
                dolfin::log (dolfin::DBG, "Printing results to file...");
                print (OUTFILE, 
                       minimizationIteration, 
                       currentFunctionalValue, 
                       alpha, 
                       backtrackingIteration, 
                       gradientNorm, 
                       relativeIncrement);
            }
            
            dolfin::end (); // "Minimization iteration %d"
        }
        
        dolfin::end (); // "Minimization loop"
        
        
        if (isConverged () == false) 
        {
            dolfin::log (dolfin::INFO, "End of Minimization loop");
            dolfin::warning ("Maximum number of iterations reached");
            dolfin::log (dolfin::INFO, "Iterations performed: %d", minimizationIteration);
            dolfin::log (dolfin::INFO, "Gradient norm = %f", gradientNorm);
            dolfin::log (dolfin::INFO, "Relative increment = %f", relativeIncrement);
            dolfin::log (dolfin::INFO, "Functional value = %f\n", currentFunctionalValue);
        }
        else
        {
            dolfin::log (dolfin::INFO, "End of Minimization loop");
            dolfin::log (dolfin::INFO, "Iterations performed: %d", minimizationIteration);
            dolfin::log (dolfin::INFO, "Gradient norm = %f", gradientNorm);
            dolfin::log (dolfin::INFO, "Relative increment = %f", relativeIncrement);
            dolfin::log (dolfin::INFO, "Functional value = %f\n", currentFunctionalValue);
        }
        
        dolfin::log (dolfin::INFO, "-----------------------------\n");
        
        if (hasOutputFile)
        {
            OUTFILE.close ();
        }
    }
     
    

    // ***************************************** //
    // ********** PRIVATE MEMBERS ************** //
    // ***************************************** //
    bool BacktrackingOptimizer:: backtrackingLoop (const double& previousFunctionalValue,
                                                   double& currentFunctionalValue, 
                                                   const double& gradientDotSearchDirection,
                                                   double& alpha,
                                                   int& backtrackingIteration,
                                                   dolfin::Function& controlVariable,
                                                   const dolfin::Function& previousControlVariable,
                                                   const dolfin::Function& searchDirection,
                                                   dcp::GenericEquationSystem& problem,
                                                   const dcp::GenericObjectiveFunctional& objectiveFunctional, 
                                                   const dcp::GenericDescentMethod::Updater& updater)
    {
        // get parameters
        double c_1 = this->parameters ["c_1"];
        double rho = this->parameters ["rho"];
        int maxBacktrackingIterations = this->parameters ["max_backtracking_iterations"];
        
        // function to check if sufficient-decrease condition is satisfied
        auto sufficientDecreaseConditionIsSatisfied = [&c_1, &gradientDotSearchDirection] 
            (const double& previousValue, const double& currentValue, const double& alpha) 
        {
            return currentValue < previousValue + c_1 * alpha * gradientDotSearchDirection;
        };
        
        dolfin::begin (dolfin::INFO, "Starting backtracking loop...");
        while (sufficientDecreaseConditionIsSatisfied (previousFunctionalValue, currentFunctionalValue, alpha) == false 
               && backtrackingIteration < maxBacktrackingIterations) 
        {
            // iteration-specific variable initialization
            ++backtrackingIteration;
            alpha = alpha * rho;
    
            dolfin::log (dolfin::INFO, "========== Backtracking iteration %d ==========", backtrackingIteration);
            dolfin::log (dolfin::INFO, "Alpha = %f", alpha);
    
            // update control variable value
            dolfin::begin (dolfin::PROGRESS, "Updating control variable...");
            controlVariable = previousControlVariable + (searchDirection * alpha);
            dolfin::end ();
    
            // update problem 
            dolfin::begin (dolfin::PROGRESS, "Updating differential problem...");
            updater (problem, controlVariable);
            dolfin::end ();
    
            // solve problem
            dolfin::begin (dolfin::PROGRESS, "Solving differential problem...");
            problem.solve ();
            dolfin::end ();
    
            // update value of the functional
            dolfin::begin (dolfin::PROGRESS, "Evaluating functional...");
            currentFunctionalValue = objectiveFunctional.evaluateFunctional ();
            dolfin::end ();
    
            dolfin::log (dolfin::INFO, "Functional value = %f", currentFunctionalValue);
    
            dolfin::log (dolfin::INFO, "");
        }
       
        dolfin::end ();
        
        if (sufficientDecreaseConditionIsSatisfied (previousFunctionalValue, currentFunctionalValue, alpha) == false) 
        {
            dolfin::log (dolfin::INFO, "");
            dolfin::warning ("Backtracking loop ended because maximum number of iterations was reached");
            dolfin::log (dolfin::INFO, "Alpha (determined with backtracking) = %f", alpha);
            return false;
        }
        else
        {
            dolfin::log (dolfin::INFO, "");
            dolfin::log (dolfin::INFO, "Backtracking loop ended. Iterations performed: %d\n", backtrackingIteration);
            dolfin::log (dolfin::INFO, "Alpha (determined with backtracking) = %f", alpha);
            return true;
        }
    }
    
    
    
    void BacktrackingOptimizer::print (std::ostream& OUTSTREAM, 
                                       const int& iteration,
                                       const double& functionalValue,
                                       const double& alpha,
                                       const int& backtrackingIterations,
                                       const double& gradientNorm,
                                       const double& relativeIncrement)
    {
        dolfin::log (dolfin::DBG, "Printing results to file...");
        std::string convergenceCriterion = this->parameters ["convergence_criterion"];

        if (convergenceCriterion == "both")
        {
            OUTSTREAM << "  "
                << std::left
                << std::setw (14)
                << iteration 
                << std::setw (21)
                << functionalValue 
                << std::setw (15)
                << alpha
                << std::setw (28)
                << backtrackingIterations
                << std::setw (18)
                << gradientNorm 
                << std::setw (23)
                << relativeIncrement
                << std::endl;
        }
        else if (convergenceCriterion == "increment")
        {
            OUTSTREAM << "  "
                << std::left
                << std::setw (14)
                << iteration 
                << std::setw (21)
                << functionalValue 
                << std::setw (15)
                << alpha
                << std::setw (28)
                << backtrackingIterations
                << std::setw (23)
                << relativeIncrement
                << std::endl;
        }
        else if (convergenceCriterion == "gradient")
        {
            OUTSTREAM << "  "
                << std::left
                << std::setw (14)
                << iteration 
                << std::setw (21)
                << functionalValue 
                << std::setw (15)
                << alpha
                << std::setw (28)
                << backtrackingIterations
                << std::setw (18)
                << gradientNorm 
                << std::endl;
        }
        else
        {
            dolfin::dolfin_error ("dcp: BacktrackingOptimizer.cpp",
                                  "print",
                                  "Unknown convergence criterion \"%s\"", 
                                  convergenceCriterion.c_str ());
        }
    }
    


    bool BacktrackingOptimizer::openOutputFile (std::ofstream& OUTFILE)
    {
        std::string outputFileName = this->parameters ["output_file_name"];
        std::string convergenceCriterion = this->parameters ["convergence_criterion"];
        
        if (!outputFileName.empty ())
        {
            OUTFILE.open (outputFileName);
            if (OUTFILE.fail ())
            {
                dolfin::dolfin_error ("dcp: BacktrackingOptimizer.cpp",
                                      "openOutputFile",
                                      "Cannot open output file \"%s\"", 
                                      outputFileName.c_str ());
            }
            
            if (convergenceCriterion == "both")
            {
                OUTFILE << "# Iteration"
                        << "     "
                        << "Functional_value" 
                        << "     "
                        << "Alpha" 
                        << "          "
                        << "Backtracking_iterations"
                        << "     "
                        << "Gradient_norm" 
                        << "     "
                        << "Relative_increment" 
                        << std::endl;
            }
            else if (convergenceCriterion == "increment")
            {
                OUTFILE << "# Iteration"
                        << "     "
                        << "Functional_value" 
                        << "     "
                        << "Alpha" 
                        << "          "
                        << "Backtracking_iterations"
                        << "     "
                        << "Relative_increment" 
                        << std::endl;
            }
            else if (convergenceCriterion == "gradient")
            {
                OUTFILE << "# Iteration"
                        << "     "
                        << "Functional_value" 
                        << "Alpha" 
                        << "          "
                        << "Backtracking_iterations"
                        << "     "
                        << "Gradient_norm" 
                        << std::endl;
            }
            else
            {
                dolfin::dolfin_error ("dcp: BacktrackingOptimizer.cpp",
                                      "openOutputFile",
                                      "Unknown convergence criterion \"%s\"", 
                                      convergenceCriterion.c_str ());
                
            }
            return true;
        }
        
        return false;
    }
}
