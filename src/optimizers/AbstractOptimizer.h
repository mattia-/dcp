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

#ifndef SRC_OPTIMIZERS_ABSTRACTOPTIMIZER_HPP_INCLUDE_GUARD
#define SRC_OPTIMIZERS_ABSTRACTOPTIMIZER_HPP_INCLUDE_GUARD

#include <dolfin/parameter/Parameters.h>
#include <dolfin/function/GenericFunction.h>
#include <dolfin/function/Function.h>
#include <dolfin/function/Expression.h>
#include <objective_functional/AbstractObjectiveFunctional.h>
#include <differential_problems/CompositeDifferentialProblem.h>
#include <functional>

namespace dcp
{
    /*! \class AbstractOptimizer AbstractOptimizer.h
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
             *  \li the dolfin function that will contain the search direction after the function exits
             *  \li the dolfin function containing the gradient
             */
            virtual void apply (dcp::CompositeDifferentialProblem& problem,
                                const dcp::AbstractObjectiveFunctional& objectiveFunctional, 
                                dolfin::Function& initialGuess,
                                const std::function 
                                <
                                    void (dcp::CompositeDifferentialProblem&, const dolfin::GenericFunction&)
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
