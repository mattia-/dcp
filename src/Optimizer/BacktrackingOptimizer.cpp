#include <Optimizer/BacktrackingOptimizer.hpp>
#include <Optimizer/dotproduct.h>
#include <dolfin/parameter/Parameters.h>
#include <dolfin/log/dolfin_log.h>
#include <dolfin/function/Function.h>
#include <dolfin/common/NoDeleter.h>
#include <dolfin/fem/assemble.h>
#include <dolfin/fem/Form.h>
#include <boost/shared_ptr.hpp>
#include <string>
#include <functional>
#include <cmath>
#include <iomanip>

namespace controlproblem
{
    BacktrackingOptimizer::BacktrackingOptimizer (const double& gradientNormTolerance,
                                                  const double& relativeIncrementTolerance,
                                                  const std::string& convergenceCriterion) :
        AbstractOptimizer (),
        dotProductComputer_ ()
    {
        dolfin::begin (dolfin::DBG, "Creating BacktrackingOptimizer object...");
        
        dolfin::log (dolfin::DBG, "Setting up parameters...");
        parameters.add ("descent_method", "backtracking_gradient_method");
        parameters.add ("gradient_norm_tolerance", gradientNormTolerance);
        parameters.add ("relative_increment_tolerance", relativeIncrementTolerance);
        parameters.add ("convergence_criterion", convergenceCriterion);
        parameters.add ("c_1", 1e-3);
        parameters.add ("alpha_0", 0.5);
        parameters.add ("rho", 0.5);
        parameters.add ("max_minimization_iterations", 100);
        parameters.add ("max_backtracking_iteration", 20);
        parameters.add ("output_file_name", "");
        
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
        double c_1                        = this->parameters ["c_1"];
        double alpha_0                    = this->parameters ["alpha_0"];
        double rho                        = this->parameters ["rho"];
        double gradientNormTolerance      = this->parameters ["gradient_norm_tolerance"];
        double relativeIncrementTolerance = this->parameters ["relative_increment_tolerance"];
        std::string convergenceCriterion  = this->parameters ["convergence_criterion"];
        int maxMinimizationIterations     = this->parameters ["max_minimization_iterations"];
        int maxBacktrackingIterations     = this->parameters ["max_backtracking_iteration"];
        std::string outputFileName        = this->parameters ["output_file_name"];

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
        boost::shared_ptr<dolfin::Form> dotProductComputer;
        
        // define output file and print header if necessary
        bool hasOutputFile = false;
        std::ofstream OUTFILE;
        if (!outputFileName.empty ())
        {
            OUTFILE.open (outputFileName);
            if (OUTFILE.fail ())
            {
                dolfin::error ("Cannot open output file \"%s\"", outputFileName.c_str ());
            }
            
            hasOutputFile = true;
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
                dolfin::error ("Unknown convergence criterion \"%s\"", convergenceCriterion.c_str ());
            }
        }
        
        
        // all linear problems in the CompositeDifferentialProblem "problem" should be reassembled every time. So we
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
        

        // ------------------------------------------------------------------------------------------------------- //
        // create correct object-type for dotProductComputer
        if (dotProductComputer_ != nullptr)
        {
            dolfin::log (dolfin::DBG, "Using protected member variable to compute dot products and norms");
            dotProductComputer = dotProductComputer_;
        }
        else
        {
            // get type of cells in mesh
            std::string meshCellType = dolfin::CellType::type2string ((controlVariable.function_space ()->mesh ())->type ().cell_type ());
            dolfin::log (dolfin::DBG, "Mesh cell type is: %s", meshCellType.c_str ());

            // get rank of control function
            int controlVariableRank = controlVariable.value_rank ();
            dolfin::log (dolfin::DBG, "Control function rank is: %d", controlVariableRank);
            
            // get control variable mesh
            boost::shared_ptr<const dolfin::Mesh> controlVariableMesh = controlVariable.function_space () -> mesh ();
 
            if (meshCellType == "interval")
            {
                if (controlVariableRank == 0)
                {
                    dolfin::log (dolfin::DBG, "Selected scalar 1D form to compute dot products and norms");
                    dotProductComputer.reset (new dotproduct::Form_scalar1D_dotProduct (controlVariableMesh));
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
                    dotProductComputer.reset (new dotproduct::Form_scalar2D_dotProduct (controlVariableMesh));
                }
                else if (controlVariableRank == 1)
                {
                    dolfin::log (dolfin::DBG, "Selected vector 2D form to compute dot products and norms");
                    dotProductComputer.reset (new dotproduct::Form_vector2D_dotProduct (controlVariableMesh));
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
                    dotProductComputer.reset (new dotproduct::Form_scalar3D_dotProduct (controlVariableMesh));
                }
                else if (controlVariableRank == 1)
                {
                    dolfin::log (dolfin::DBG, "Selected vector 3D form to compute dot products and norms");
                    dotProductComputer.reset (new dotproduct::Form_vector3D_dotProduct (controlVariableMesh));
                }
                else
                {
                    dolfin::error ("No form to compute dot products and norms for mesh cell type \"%s\" and control function rank %d", 
                                   meshCellType.c_str (),
                                   controlVariableRank);
                }
            }
        }
        // end of if statement to correctly set dotProductComputer
        // ------------------------------------------------------------------------------------------------------- //
        

        // ------------------------------------------------------------------------------------------------------- //
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
            dolfin::error ("Unknown convergence criterion \"%s\"", convergenceCriterion.c_str ());
        }
        
        
        // function to check if sufficient-decrease condition is satisfied
        auto sufficientDecreaseConditionIsSatisfied = [&c_1, &gradientDotSearchDirection] 
            (const double& previousValue, const double& currentValue, const double& alpha) 
        {
            return currentValue < previousValue + c_1 * alpha * gradientDotSearchDirection;
        };
        // end of minimization loop functions setting
        // ------------------------------------------------------------------------------------------------------- //

        
        // update and solve the problem for the first time
        dolfin::begin ("Minimization loop initialization...");
        
        dolfin::log (dolfin::PROGRESS, "Updating differential problem using control initial guess...");
        updater (problem, controlVariable);

        dolfin::log (dolfin::PROGRESS, "Solving differential problem...");
        problem.solve ();
        
        
        // initialize loop variables
        dotProductComputer -> set_coefficient (0, dolfin::reference_to_no_delete_pointer (objectiveFunctional.gradient ()));
        dotProductComputer -> set_coefficient (1, dolfin::reference_to_no_delete_pointer (objectiveFunctional.gradient ()));
        gradientNorm = sqrt (dolfin::assemble (*dotProductComputer));
        
        relativeIncrement = relativeIncrementTolerance + 1; // just an initialization to be sure that the first iteration 
                                                            // of the minimization loop is performed
        currentFunctionalValue = objectiveFunctional.evaluateFunctional ();


        dolfin::log (dolfin::DBG, "======== MINIMIZATION LOOP ========");
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
        dolfin::log (dolfin::INFO, "======== LOOP VARIABLES INITIAL VALUE ========");
        dolfin::log (dolfin::INFO, "gradient norm = %f", gradientNorm);
        dolfin::log (dolfin::INFO, "functional value = %f\n", currentFunctionalValue);
            
        
        dolfin::begin ("Starting minimization loop...");
        int minimizationIteration = 0;
        
        if (hasOutputFile)
        {
            print (OUTFILE, minimizationIteration, currentFunctionalValue, 0, 0, gradientNorm, relativeIncrement);
        }
        
        while (isConverged () == false && minimizationIteration < maxMinimizationIterations)
        {
            minimizationIteration++;
            
            dolfin::log (dolfin::INFO, "=========================");
            dolfin::log (dolfin::INFO, "Minimization iteration %d", minimizationIteration);
            dolfin::begin ("=========================");
            
            // iteration-specific variable initialization
            alpha = alpha_0;
            int backtrackingIteration = 0;
            previousFunctionalValue = currentFunctionalValue;
            functionalGradient = objectiveFunctional.gradient ();
            searchDirectionComputer (searchDirection, functionalGradient);
            
            // compute dot product between gradient and search direction
            dotProductComputer -> set_coefficient (0, dolfin::reference_to_no_delete_pointer (objectiveFunctional.gradient ()));
            dotProductComputer -> set_coefficient (1, dolfin::reference_to_no_delete_pointer (searchDirection));
            gradientDotSearchDirection = dolfin::assemble (*dotProductComputer);
            
            // solution of problem with alpha_0
            dolfin::log (dolfin::INFO, "Alpha = %f", alpha);

            // update control variable value
            dolfin::log (dolfin::PROGRESS, "Updating control variable...");
            dolfin::Function previousControlVariable (controlVariable);
            controlVariable = previousControlVariable + (searchDirection * alpha);

            // update problem 
            dolfin::log (dolfin::PROGRESS, "Updating differential problem...");
            updater (problem, controlVariable);

            // solve problem
            dolfin::log (dolfin::PROGRESS, "Solving differential problem...");
            problem.solve ();

            // update value of the functional
            dolfin::log (dolfin::PROGRESS, "Evaluating functional...");
            currentFunctionalValue = objectiveFunctional.evaluateFunctional ();
            
            dolfin::log (dolfin::INFO, "Functional value = %f\n", currentFunctionalValue);
            
            
            dolfin::begin (dolfin::INFO, "Starting backtracking loop...");
            
            // backtacking loop
            while (sufficientDecreaseConditionIsSatisfied (previousFunctionalValue, currentFunctionalValue, alpha) == false 
                   && backtrackingIteration < maxBacktrackingIterations) 
            {
                // iteration-specific variable initialization
                ++backtrackingIteration;
                alpha = alpha * rho;
                
                dolfin::log (dolfin::INFO, "========== Backtracking iteration %d ==========", backtrackingIteration);
                dolfin::log (dolfin::INFO, "Alpha = %f", alpha);
                
                // update control variable value
                dolfin::log (dolfin::PROGRESS, "Updating control variable...");
                controlVariable = previousControlVariable + (searchDirection * alpha);
                
                // update problem 
                dolfin::log (dolfin::PROGRESS, "Updating differential problem...");
                updater (problem, controlVariable);
                
                // solve problem
                dolfin::log (dolfin::PROGRESS, "Solving differential problem...");
                problem.solve ();
                
                // update value of the functional
                dolfin::log (dolfin::PROGRESS, "Evaluating functional...");
                currentFunctionalValue = objectiveFunctional.evaluateFunctional ();
                
                dolfin::log (dolfin::INFO, "Functional value = %f", currentFunctionalValue);
                
                dolfin::log (dolfin::INFO, "");
            }
            
            dolfin::end ();
             
            if (sufficientDecreaseConditionIsSatisfied (previousFunctionalValue, currentFunctionalValue, alpha) == false) 
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
                
            
            // update gradient norm if necessary for convergence check
            if (convergenceCriterion == "gradient" || convergenceCriterion == "both")
            {
                dolfin::log (dolfin::DBG, "Computing gradient norm...");
                dotProductComputer -> set_coefficient (0, dolfin::reference_to_no_delete_pointer (objectiveFunctional.gradient ()));
                dotProductComputer -> set_coefficient (1, dolfin::reference_to_no_delete_pointer (objectiveFunctional.gradient ()));
                gradientNorm = sqrt (dolfin::assemble (*dotProductComputer));
                dolfin::log (dolfin::INFO, "Gradient norm = %f", gradientNorm);
            }
            

            // update relative increment if necessary for convergence check
            if (convergenceCriterion == "increment" || convergenceCriterion == "both")
            {
                dolfin::log (dolfin::DBG, "Computing relative increment...");
                dotProductComputer -> set_coefficient (0, dolfin::reference_to_no_delete_pointer (previousControlVariable));
                dotProductComputer -> set_coefficient (1, dolfin::reference_to_no_delete_pointer (previousControlVariable));
                previousControlVariableNorm = sqrt (dolfin::assemble (*dotProductComputer));

                controlVariableIncrement = controlVariable - previousControlVariable;
                dotProductComputer -> set_coefficient (0, dolfin::reference_to_no_delete_pointer (controlVariableIncrement));
                dotProductComputer -> set_coefficient (1, dolfin::reference_to_no_delete_pointer (controlVariableIncrement));
                incrementNorm = sqrt (dolfin::assemble (*dotProductComputer));

                // compute relative increment. We add DOLFIN_EPS at the denominator in case previousControlVariableNorm 
                // is zero (like at the first iteration)
                relativeIncrement = incrementNorm / (previousControlVariableNorm + DOLFIN_EPS); 
                dolfin::log (dolfin::INFO, "Relative increment = %f", relativeIncrement);
            }
                
            dolfin::log (dolfin::INFO, "Functional value = %f\n\n", currentFunctionalValue);
            
            if (hasOutputFile)
            {
                print (OUTFILE, 
                       minimizationIteration, 
                       currentFunctionalValue, 
                       alpha, 
                       backtrackingIteration, 
                       gradientNorm, 
                       relativeIncrement);
            }
        }
        
        
        if (isConverged () == false) 
        {
            dolfin::log (dolfin::INFO, "");
            dolfin::warning ("Minimization loop ended because maximum number of iterations was reached");
            dolfin::log (dolfin::INFO, "Gradient norm = %f", gradientNorm);
            dolfin::log (dolfin::INFO, "Relative increment = %f", relativeIncrement);
            dolfin::log (dolfin::INFO, "Functional value = %f\n\n", currentFunctionalValue);
        }
        else
        {
            dolfin::log (dolfin::INFO, "Minimization loop ended. Iterations performed: %d\n", minimizationIteration);
            dolfin::log (dolfin::INFO, "Gradient norm = %f", gradientNorm);
            dolfin::log (dolfin::INFO, "Relative increment = %f", relativeIncrement);
            dolfin::log (dolfin::INFO, "Functional value = %f\n\n", currentFunctionalValue);
        }
        
        dolfin::end ();
        
        if (hasOutputFile)
        {
            OUTFILE.close ();
        }
    }
    


    void BacktrackingOptimizer::gradientSearchDirection (dolfin::Function& searchDirection, const dolfin::Function& gradient)
    {
        searchDirection = gradient * (-1);
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
            dolfin::error ("Unknown convergence criterion \"%s\"", convergenceCriterion.c_str ());
        }
    }
}
