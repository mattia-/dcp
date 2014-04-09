#include <Optimizer/BacktrackingOptimizer.hpp>
#include <dolfin/parameter/Parameters.h>
#include <dolfin/log/dolfin_log.h>
#include <dolfin/function/Function.h>
#include <dolfin/common/NoDeleter.h>
#include <dolfin/fem/assemble.h>
#include <string>
#include <functional>
#include <cmath>
#include <Optimizer/AuxiliaryForms.h>

namespace controlproblem
{
    BacktrackingOptimizer::BacktrackingOptimizer (const double& gradientNormTolerance,
                                                  const double& incrementNormTolerance,
                                                  const std::string& convergenceCriterion) :
        AbstractOptimizer ()
    {
        dolfin::begin (dolfin::DBG, "Creating BacktrackingOptimizer object...");
        
        dolfin::log (dolfin::DBG, "Setting up parameters...");
        parameters.add ("descent_method", "backtracking_gradient_method");
        parameters.add ("gradient_norm_tolerance", gradientNormTolerance);
        parameters.add ("increment_norm_tolerance", incrementNormTolerance);
        parameters.add ("convergence_criterion", convergenceCriterion);
        parameters.add ("c_1", 1e-3);
        parameters.add ("alpha_0", 0.5);
        parameters.add ("rho", 0.5);
        parameters.add ("max_minimization_iterations", 100);
        parameters.add ("max_alpha_iterations", 10);
        
        dolfin::log (dolfin::DBG, "BacktrackingOptimizer object created");
        
        dolfin::end ();
    }
    
    

    void BacktrackingOptimizer::apply (controlproblem::CompositeDifferentialProblem& problem,
                                       const controlproblem::AbstractObjectiveFunctional& objectiveFunctional, 
                                       dolfin::Function& initialGuess,
                                       const std::function 
                                       <
                                           void (controlproblem::CompositeDifferentialProblem&, const dolfin::GenericFunction&)
                                       >& updater)
    {
        // get parameters values
        double c_1                       = this->parameters ["c_1"];
        double alpha_0                   = this->parameters ["alpha_0"];
        double rho                       = this->parameters ["rho"];
        double gradientNormTolerance     = this->parameters ["gradient_norm_tolerance"];
        double incrementNormTolerance    = this->parameters ["increment_norm_tolerance"];
        std::string convergenceCriterion = this->parameters ["convergence_criterion"];
        int maxMinimizationIterations    = this->parameters ["max_minimization_iterations"];
        int maxBacktrackingIterations           = this->parameters ["max_alpha_iterations"];

        // define loop variables
        double alpha;
        double gradientNorm;
        double incrementNorm;
        double functionalValueOld;
        double functionalValue;
        dolfin::Function controlVariable (initialGuess);
        dolfin::Function functionalGradient (controlVariable.function_space ());

        // define variable to compute the norm of the gradient of the functional
        AuxiliaryForms::Form_NormComputer normComputer (objectiveFunctional.mesh ());
        normComputer.set_coefficient ("u", dolfin::reference_to_no_delete_pointer (objectiveFunctional.gradient ()));
        
        // function to check convergence
        std::function <bool ()> isConverged;
        
        // set value of isConverged () according to convergence criterion
        if (convergenceCriterion == "gradient_norm")
        {
            isConverged = [&] () {return gradientNorm < gradientNormTolerance;};
        }
        else if (convergenceCriterion == "increment_norm")
        {
            isConverged = [&] () {return incrementNorm < incrementNormTolerance;};
        }
        else if (convergenceCriterion == "both")
        {
            isConverged = [&] () {return (incrementNorm < incrementNormTolerance) || (gradientNorm < gradientNormTolerance);};
        }
        else
        {
            dolfin::error ("Unknown convergence criterion \"%s\"", convergenceCriterion.c_str ());
        }
        
        // function to check if sufficient-decrease condition is satisfied
        auto sufficientDecreaseConditionIsSatisfied = [&c_1, &gradientNorm] 
            (const double& oldValue, const double& newValue, const double& alpha) 
        {
            return newValue < oldValue + c_1 * alpha * gradientNorm;
        };

        // update and solve the problem for the first time
        updater (problem, controlVariable);
        problem.solve ();
        
        // initialize loop variables
        alpha = alpha_0;
        gradientNorm = dolfin::assemble (normComputer);
        incrementNorm = incrementNormTolerance + 1; // just an initialization to be sure that the first iteration 
                                                    // of the minimization loop is performed
        functionalValue = 0.0;
        functionalValueOld = functionalValue;

        dolfin::log (dolfin::DBG, "======== MINIMIZATION LOOP ========");
        dolfin::log (dolfin::DBG, "Parameters value are:");
        dolfin::log (dolfin::DBG, "c_1 = %f", c_1);
        dolfin::log (dolfin::DBG, "alpha_0 = %f", alpha_0);
        dolfin::log (dolfin::DBG, "rho = %f", rho);
        dolfin::log (dolfin::DBG, "convergence check tolerance on gradient norm = %f", gradientNormTolerance);
        dolfin::log (dolfin::DBG, "convergence check tolerance on increment norm = %f", incrementNormTolerance);
        dolfin::log (dolfin::DBG, "convergence criterion = %f", convergenceCriterion.c_str ());
        dolfin::log (dolfin::DBG, "maximum minimization iterations = %f", maxMinimizationIterations);
        dolfin::log (dolfin::DBG, "maximum backtracking iterations = %f", maxBacktrackingIterations);
        dolfin::log (dolfin::DBG, "gradient norm = %f", gradientNorm);
        dolfin::log (dolfin::DBG, "increment norm = %f", incrementNorm);
            
        dolfin::begin ("Starting minimization loop...");
        int minimizationIteration = 0;
        while (isConverged () == false && minimizationIteration < maxMinimizationIterations)
        {
            minimizationIteration++;
            
            dolfin::begin ("Iteration %d", minimizationIteration);
            
            // iteration-specific variable initialization
            alpha = alpha_0;
            int backtrackingIteration = 0;
            functionalValueOld = functionalValue;
            
            dolfin::begin (dolfin::PROGRESS, "Starting backtracking loop...");
            
            // note that any call to sufficientDecreaseConditionIsSatisfied () at this point will use the gradient norm
            // at the previous iteration, i.e. the gradient norm at alpha = 0
            while (sufficientDecreaseConditionIsSatisfied (functionalValueOld, functionalValue, alpha) == false 
                   && backtrackingIteration < maxBacktrackingIterations) 
            {
                backtrackingIteration++;
                if (backtrackingIteration > 1)
                {
                    alpha = alpha * rho;
                }
                
                dolfin::log (dolfin::PROGRESS, "Iteration %d", minimizationIteration);
                dolfin::log (dolfin::PROGRESS, "Alpha = %f", alpha);
                
                // update control variable value
                dolfin::log (dolfin::PROGRESS, "Updating control variable...");
                functionalGradient = objectiveFunctional.gradient ();
                controlVariable = controlVariable + functionalGradient * alpha;
                
                // update problem 
                dolfin::log (dolfin::PROGRESS, "Updating differential problem...");
                updater (problem, controlVariable);
                 
                // solve problem
                dolfin::log (dolfin::PROGRESS, "Updating differential problem...");
                problem.solve ();
                
                // update value of the functional
                dolfin::log (dolfin::PROGRESS, "Evaluating functional...");
                functionalValue = objectiveFunctional.evaluateFunctional ();
            }
            
            dolfin::end ();
             
            // update gradient norm
            gradientNorm = assemble (normComputer);

            // update increment norm if it is used for the convergence check
            incrementNorm = fabs (functionalValue - functionalValueOld); 
            
            dolfin::log (dolfin::INFO, "Alpha (determined with backtracking) = %f", alpha);
            dolfin::log (dolfin::INFO, "Gradient norm = %f", gradientNorm);
            dolfin::log (dolfin::INFO, "Increment norm = %f", incrementNorm);
        }
        dolfin::end ();
    }
}
