#ifndef SRC_DESCENTMETHOD_ABSTRACTDESCENTMETHOD_HPP_INCLUDE_GUARD
#define SRC_DESCENTMETHOD_ABSTRACTDESCENTMETHOD_HPP_INCLUDE_GUARD

#include <dolfin/parameter/Parameters.h>
#include <dolfin/function/Function.h>
#include <dolfin/function/Expression.h>
#include <ObjectiveFunctional/AbstractObjectiveFunctional.hpp>
#include <DifferentialProblem/CompositeDifferentialProblem.hpp>

namespace controlproblem
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
             *  \param point the starting point for the minimization algorithm
             */
            virtual void apply (const controlproblem::CompositeDifferentialProblem& problem,
                                const controlproblem::AbstractObjectiveFunctional& objectiveFunctional, 
                                dolfin::Function& point) = 0;
            

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
