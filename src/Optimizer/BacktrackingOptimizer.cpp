#include <Optimizer/BacktrackingOptimizer.hpp>
#include <dolfin/parameter/Parameters.h>
#include <dolfin/log/dolfin_log.h>
#include <dolfin/function/Function.h>
#include <dolfin/common/NoDeleter.h>
#include <dolfin/fem/assemble.h>
#include <dolfin/fem/Form.h>
#include <string>
#include <functional>
#include <cmath>
#include <Optimizer/dotproduct.h>
#include <boost/shared_ptr.hpp>

namespace controlproblem
{
    BacktrackingOptimizer::BacktrackingOptimizer (const double& gradientNormTolerance,
                                                  const double& incrementNormTolerance,
                                                  const std::string& convergenceCriterion) :
        AbstractOptimizer (),
        dotProductComputer_ ()
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
        parameters.add ("max_backtracking_iteration", 100);
        
        dolfin::log (dolfin::DBG, "BacktrackingOptimizer object created");
        
        dolfin::end ();
    }
    


    void BacktrackingOptimizer::setDotProductComputer (const boost::shared_ptr<dolfin::Form> dotProductComputer)
    {
        dotProductComputer_ = dotProductComputer;
    }



    void BacktrackingOptimizer::resetDotProductComputer ()
    {
        dotProductComputer_.reset ();
    }

    

    void BacktrackingOptimizer::apply (controlproblem::CompositeDifferentialProblem& problem,
                                       const controlproblem::AbstractObjectiveFunctional& objectiveFunctional, 
                                       dolfin::Function& initialGuess,
                                       const std::function 
                                       <
                                           void (controlproblem::CompositeDifferentialProblem&, const dolfin::GenericFunction&)
                                       >& updater,
                                       const std::function
                                       <
                                           void (dolfin::Function&, const dolfin::Function&)
                                       >& searchDirectionComputer)
    {
        // get parameters values
        double c_1                       = this->parameters ["c_1"];
        double alpha_0                   = this->parameters ["alpha_0"];
        double rho                       = this->parameters ["rho"];
        double gradientNormTolerance     = this->parameters ["gradient_norm_tolerance"];
        double incrementNormTolerance    = this->parameters ["increment_norm_tolerance"];
        std::string convergenceCriterion = this->parameters ["convergence_criterion"];
        int maxMinimizationIterations    = this->parameters ["max_minimization_iterations"];
        int maxBacktrackingIterations    = this->parameters ["max_backtracking_iteration"];

        // define loop variables
        double alpha;
        double gradientNorm;
        double gradientDotSearchDirection;
        double incrementNorm;
        double functionalValueOld;
        double functionalValue;
        dolfin::Function& controlVariable = initialGuess;
        dolfin::Function functionalGradient (controlVariable.function_space ());
        dolfin::Function searchDirection (controlVariable.function_space ());
        boost::shared_ptr<dolfin::Form> dotProductComputer;
        boost::shared_ptr<dolfin::Form> normComputer;
        

        // create correct objects for dotProductComputer and normComputer
        if (dotProductComputer != nullptr)
        {
            dolfin::log (dolfin::DBG, "Using input argument form to compute dot products and norms");
            dotProductComputer = dotProductComputer_;
            normComputer = dotProductComputer_;
        }
        else
        {
            // get type of cells in mesh
            std::string meshCellType = dolfin::CellType::type2string ((controlVariable.function_space ()->mesh ())->type ().cell_type ());
            dolfin::log (dolfin::DBG, "Mesh cell type is: %s", meshCellType.c_str ());

            // get rank of control function
            int controlVariableRank = controlVariable.value_rank ();
            dolfin::log (dolfin::DBG, "Control function rank is: %d", controlVariableRank);

            if (meshCellType == "interval")
            {
                if (controlVariableRank == 0)
                {
                    dolfin::log (dolfin::DBG, "Selected scalar 1D form to compute dot products and norms");
                    dotProductComputer.reset (new dotproduct::Form_scalar1D_dotProduct (objectiveFunctional.mesh ()));
                    normComputer.reset (new dotproduct::Form_scalar1D_dotProduct (objectiveFunctional.mesh ()));
                }
                else
                {
                    dolfin::error ("No form to compute dot products and norms for mesh cell type \"%s\" and control function rank %d", 
                                   meshCellType.c_str (),
                                   controlVariableRank);
                }
            }

            if (meshCellType == "triangle")
            {
                if (controlVariableRank == 0)
                {
                    dolfin::log (dolfin::DBG, "Selected scalar 2D form to compute dot products and norms");
                    dotProductComputer.reset (new dotproduct::Form_scalar2D_dotProduct (objectiveFunctional.mesh ()));
                    normComputer.reset (new dotproduct::Form_scalar2D_dotProduct (objectiveFunctional.mesh ()));
                }
                if (controlVariableRank == 1)
                {
                    dolfin::log (dolfin::DBG, "Selected vectorial 2D form to compute dot products and norms");
                    dotProductComputer.reset (new dotproduct::Form_vectorial2D_dotProduct (objectiveFunctional.mesh ()));
                    normComputer.reset (new dotproduct::Form_vectorial2D_dotProduct (objectiveFunctional.mesh ()));
                }
                else
                {
                    dolfin::error ("No form to compute dot products and norms for mesh cell type \"%s\" and control function rank %d", 
                                   meshCellType.c_str (),
                                   controlVariableRank);
                }
            }

            if (meshCellType == "tetrahedron")
            {
                if (controlVariableRank == 0)
                {
                    dolfin::log (dolfin::DBG, "Selected scalar 3D form to compute dot products and norms");
                    dotProductComputer.reset (new dotproduct::Form_scalar3D_dotProduct (objectiveFunctional.mesh ()));
                    normComputer.reset (new dotproduct::Form_scalar3D_dotProduct (objectiveFunctional.mesh ()));
                }
                if (controlVariableRank == 1)
                {
                    dolfin::log (dolfin::DBG, "Selected vectorial 3D form to compute dot products and norms");
                    dotProductComputer.reset (new dotproduct::Form_vectorial3D_dotProduct (objectiveFunctional.mesh ()));
                    normComputer.reset (new dotproduct::Form_vectorial3D_dotProduct (objectiveFunctional.mesh ()));
                }
                else
                {
                    dolfin::error ("No form to compute dot products and norms for mesh cell type \"%s\" and control function rank %d", 
                                   meshCellType.c_str (),
                                   controlVariableRank);
                }
            }
        }
        
        
        // set coefficients for dot product and norm computers
        normComputer -> set_coefficient (0, dolfin::reference_to_no_delete_pointer (objectiveFunctional.gradient ()));
        normComputer -> set_coefficient (1, dolfin::reference_to_no_delete_pointer (objectiveFunctional.gradient ()));
        
        dotProductComputer -> set_coefficient (0, dolfin::reference_to_no_delete_pointer (objectiveFunctional.gradient ()));
        dotProductComputer -> set_coefficient (1, dolfin::reference_to_no_delete_pointer (searchDirection));
        

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
        auto sufficientDecreaseConditionIsSatisfied = [&c_1, &gradientDotSearchDirection] 
            (const double& oldValue, const double& newValue, const double& alpha) 
        {
            return newValue < oldValue - c_1 * alpha * gradientDotSearchDirection;
        };

        
        // update and solve the problem for the first time
        dolfin::begin ("Minimization loop initialization...");
        
        dolfin::log (dolfin::PROGRESS, "Updating differential problem using control initial guess...");
        updater (problem, controlVariable);

        dolfin::log (dolfin::PROGRESS, "Solving differential problem...");
        problem.solve ();
        
        
        // initialize loop variables
        alpha = alpha_0;
        gradientNorm = dolfin::assemble (*normComputer);
        incrementNorm = incrementNormTolerance + 1; // just an initialization to be sure that the first iteration 
                                                    // of the minimization loop is performed
        functionalValue = objectiveFunctional.evaluateFunctional ();


        dolfin::log (dolfin::DBG, "======== MINIMIZATION LOOP ========");
        dolfin::log (dolfin::DBG, "Parameters value are:");
        dolfin::log (dolfin::DBG, "c_1 = %f", c_1);
        dolfin::log (dolfin::DBG, "alpha_0 = %f", alpha_0);
        dolfin::log (dolfin::DBG, "rho = %f", rho);
        dolfin::log (dolfin::DBG, "gradient norm tolerance = %f", gradientNormTolerance);
        dolfin::log (dolfin::DBG, "increment norm tolerance = %f", incrementNormTolerance);
        dolfin::log (dolfin::DBG, "convergence criterion = %s", convergenceCriterion.c_str ());
        dolfin::log (dolfin::DBG, "maximum minimization iterations = %d", maxMinimizationIterations);
        dolfin::log (dolfin::DBG, "maximum backtracking iterations = %d", maxBacktrackingIterations);
        dolfin::log (dolfin::DBG, "gradient norm = %f", gradientNorm);
        dolfin::log (dolfin::DBG, "increment norm = %f", incrementNorm);
        dolfin::log (dolfin::DBG, "functional value = %f", functionalValue);
            
        
        dolfin::begin ("Starting minimization loop...");
        int minimizationIteration = 0;
        while (isConverged () == false && minimizationIteration < maxMinimizationIterations)
        {
            minimizationIteration++;
            
            dolfin::log (dolfin::INFO, "=========================");
            dolfin::log (dolfin::INFO, "Minimization iteration %d", minimizationIteration);
            dolfin::begin ("=========================");
            
            // iteration-specific variable initialization
            alpha = alpha_0;
            int backtrackingIteration = 0;
            functionalValueOld = functionalValue;
            functionalGradient = objectiveFunctional.gradient ();
            gradientDotSearchDirection = dolfin::assemble (*dotProductComputer);
            
            // solution of problem with alpha_0
            dolfin::log (dolfin::INFO, "Alpha = %f", alpha);

            // update control variable value
            dolfin::log (dolfin::PROGRESS, "Updating control variable...");
            dolfin::Function aux (controlVariable);
            searchDirectionComputer (searchDirection, functionalGradient);
            controlVariable = aux + (searchDirection * alpha);

            // update problem 
            dolfin::log (dolfin::PROGRESS, "Updating differential problem...");
            updater (problem, controlVariable);

            // solve problem
            dolfin::log (dolfin::PROGRESS, "Solving differential problem...");
            problem.solve ();

            // update value of the functional
            dolfin::log (dolfin::PROGRESS, "Evaluating functional...");
            functionalValue = objectiveFunctional.evaluateFunctional ();
            
            dolfin::log (dolfin::INFO, "Functional value = %f", functionalValue);
            
            
            dolfin::begin (dolfin::INFO, "Starting backtracking loop...");
            
            // note that any call to sufficientDecreaseConditionIsSatisfied () at this point will use the gradient norm
            // at the previous iteration, i.e. the gradient norm at alpha = 0
            while (sufficientDecreaseConditionIsSatisfied (functionalValueOld, functionalValue, alpha) == false 
                   && backtrackingIteration < maxBacktrackingIterations) 
            {
                // iteration-specific variable initialization
                ++backtrackingIteration;
                alpha = alpha * rho;
                
                dolfin::log (dolfin::INFO, "========== Backtracking iteration %d ==========", backtrackingIteration);
                dolfin::log (dolfin::INFO, "Alpha = %f", alpha);
                
                // update control variable value
                dolfin::log (dolfin::PROGRESS, "Updating control variable...");
                searchDirectionComputer (searchDirection, functionalGradient);
                controlVariable = aux + (searchDirection * alpha);
                
                // update problem 
                dolfin::log (dolfin::PROGRESS, "Updating differential problem...");
                updater (problem, controlVariable);
                
                // solve problem
                dolfin::log (dolfin::PROGRESS, "Solving differential problem...");
                problem.solve ();
                
                // update value of the functional
                dolfin::log (dolfin::PROGRESS, "Evaluating functional...");
                functionalValue = objectiveFunctional.evaluateFunctional ();
                
                dolfin::log (dolfin::INFO, "Functional value = %f", functionalValue);
                
                dolfin::log (dolfin::INFO, "");
            }
            
            dolfin::end ();
             
            if (sufficientDecreaseConditionIsSatisfied (functionalValueOld, functionalValue, alpha) == false 
                && backtrackingIteration == maxBacktrackingIterations) 
            {
                dolfin::log (dolfin::INFO, "");
                dolfin::warning ("Backtracking loop ended because maximum number of iterations was reached");
                dolfin::log (dolfin::INFO, "Alpha (determined with backtracking) = %f", alpha);
            }
            else
            {
                dolfin::log (dolfin::INFO, "");
                dolfin::log (dolfin::INFO, "Backtracking loop ended. Iterations performed: %d\n", backtrackingIteration);
                dolfin::log (dolfin::INFO, "Alpha (determined with backtracking) = %f", alpha);
            }
                
            
            // update gradient nom
            gradientNorm = dolfin::assemble (*normComputer);

            // update increment norm if it is used for the convergence check
            incrementNorm = fabs (functionalValue - functionalValueOld); 
            
            dolfin::log (dolfin::INFO, "Gradient norm = %f", gradientNorm);
            dolfin::log (dolfin::INFO, "Increment norm = %f", incrementNorm);
            dolfin::log (dolfin::INFO, "Functional value = %f\n\n", 
                         backtrackingIteration > 0 ? functionalValue : objectiveFunctional.evaluateFunctional ());
        }
        
        dolfin::end ();
    }
    


    void BacktrackingOptimizer::gradientSearchDirection (dolfin::Function& searchDirection, const dolfin::Function& gradient)
    {
        searchDirection = gradient * (-1);
    }
}
