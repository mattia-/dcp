#ifndef SRC_OPTIMIZER_ABSTRACTOPTIMIZER_HPP_INCLUDE_GUARD
#define SRC_OPTIMIZER_ABSTRACTOPTIMIZER_HPP_INCLUDE_GUARD

#include <dolfin/parameter/Parameters.h>
#include <dolfin/function/GenericFunction.h>
#include <dolfin/function/Function.h>
#include <dolfin/function/Expression.h>
#include <ObjectiveFunctional/AbstractObjectiveFunctional.hpp>
#include <DifferentialProblem/CompositeDifferentialProblem.hpp>
#include <functional>

namespace DCP
{
    /*! \class AbstractOptimizer AbstractOptimizer.hpp
     *  \brief Abstract base class for descent methods.
     * 
     *  This class defines the base interface for all descent methods.
     *  It provides a \c apply() method to perform the optimization of the
     *  problem passed as input argument and an empty parameters set that 
     *  can be populated by derived classes to store concrete methods' settings
     */
    
    class AbstractOptimizer
    {
        // ---------------------------------------------------------------------------------------------//
        
        public:
            /************************* CONSTRUCTORS ********************/
            //! Default constructor
            AbstractOptimizer ();
            

            /************************* DESTRUCTOR ********************/
            //! Destructor
            /*! Default destructor, since members of the class are trivially 
             * destructible.
             */
            virtual ~AbstractOptimizer () {};
            
            
            /********************** METHODS ***********************/
            //! Perform optimization on the input problem
            /*! 
             *  Input arguments are:
             *  \param problem the composite differential problem that represents the primal/adjoint system
             *  \param objectiveFunctional the objective functional to be minimized
             *  \param initialGuess the starting point for the minimization algorithm. At the end of the function, it
             *  will containt the final value of the control variable
             *  \param updater callable object to update the control parameter value. It can be either be a function 
             *  pointer, a function object or a lambda expression. Its input argument are:
             *  \li the composite differential problem to update
             *  \li the new value of the control function
             *  
             *  \param searchDirectionComputer callable object to compute the search direction. It can either be a 
             *  function pointer, a function object or a lambda expression. In general the search direction is computed
             *  as: 
             *  \f[
             *      \mathbf{d}_k = -B_k\,\nabla J_k
             *  \f]
             *  The default value is the member function \c gradientSearchDirection(), that basically uses the above 
             *  formula with \f$ B_k = I \f$.
             *  The input arguments for \c searchDirectionComputer are:
             *  \li the function that will contain the search direction after the function is ended
             *  \li the function containing the gradient
             */
            virtual void apply (DCP::CompositeDifferentialProblem& problem,
                                const DCP::AbstractObjectiveFunctional& objectiveFunctional, 
                                dolfin::Function& initialGuess,
                                const std::function 
                                <
                                    void (DCP::CompositeDifferentialProblem&, const dolfin::GenericFunction&)
                                >& updater,
                                const std::function
                                <
                                    void (dolfin::Function&, const dolfin::Function&)
                                >& searchDirectionComputer) = 0;
            

            /********************** VARIABLES ***********************/
            //! the problem parameters
            dolfin::Parameters parameters;
            
            // ---------------------------------------------------------------------------------------------//

        protected:

            // ---------------------------------------------------------------------------------------------//

        private:
            
    };
}

#endif
