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

#ifndef SRC_OPTIMIZERS_ABSTRACTDESCENTMETHOD_H_INCLUDE_GUARD
#define SRC_OPTIMIZERS_ABSTRACTDESCENTMETHOD_H_INCLUDE_GUARD

#include <dolfin/parameter/Parameters.h>
#include <dolfin/function/GenericFunction.h>
#include <dolfin/function/Function.h>
#include <dolfin/function/Expression.h>
#include <dcp/objective_functional/AbstractObjectiveFunctional.h>
#include <dcp/differential_problems/AbstractEquationSystem.h>
#include <dcp/optimizers/GradientSearchDirection.h>
#include <dcp/utils/DotProduct.h>
#include <functional>

namespace dcp
{
    /*! \class AbstractDescentMethod AbstractDescentMethod.h
     *  \brief Abstract base class for descent methods.
     * 
     *  This class defines the base interface for all descent methods.
     *  It provides a \c apply() method to perform the optimization of the
     *  problem passed as input argument and an empty parameters set that 
     *  can be populated by derived classes to store concrete methods' settings
     */
    
    class AbstractDescentMethod
    {
        // ---------------------------------------------------------------------------------------------//
        
        public:
            /************************* TYPEDEFS ************************/
            typedef std::function<void (dcp::AbstractEquationSystem&, const dolfin::GenericFunction&)> Updater;
            typedef std::function<void (dolfin::Function&, const dolfin::Function&)> SearchDirectionComputer;
            
            /************************* CONSTRUCTORS ********************/
            //! Default constructor
            AbstractDescentMethod ();
            

            /************************* DESTRUCTOR ********************/
            //! Destructor
            /*! Default destructor, since members of the class are trivially 
             *  destructible.
             *  The protected member will all be initialized with default parameters (see specific members documentation
             *  for details).
             */
            virtual ~AbstractDescentMethod () {};
            
            
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
             *  Functors for the most common types of update are provided: see \c dcp::DirichletControlUpdater,
             *  \c dcp::DistributedControlUpdater and \c dcp::NeumannControlUpdater.
             *  
             */
            virtual void apply (dcp::AbstractEquationSystem& problem,
                                const dcp::AbstractObjectiveFunctional& objectiveFunctional, 
                                dolfin::Function& initialGuess,
                                const dcp::AbstractDescentMethod::Updater& updater) = 0;
            
            //! Set the dot product to be used
            /*!
             *  \param dotProductForm the form to be used when computing the dot product 
             */
            virtual void setDotProduct (const dolfin::Form& dotProductForm);
            
            //! Set the way the search direction is computed on every loop iteration during minimization
            /*!
             *  \param searchDirectionComputer the object to be used to compute the search direction
             */
            virtual void setSearchDirection (const SearchDirectionComputer& searchDirectionComputer);

            /********************** VARIABLES ***********************/
            //! the problem parameters
            dolfin::Parameters parameters;
            
            // ---------------------------------------------------------------------------------------------//

        protected:
            //! The form that will be used to compute the dot product between the gradient and the search direction. 
            /*! 
             *  The default value is on object of type \c dcp::DotProduct default-constructed, which will try to
             *  determine the right form to use by checking the geometrical dimensions of the input objects. 
             *  However, sometimes it may be useful to have a user-defined object to compute the dot product.
             *  o do so, use the function \c setDotProduct
             */
            dcp::DotProduct dotProduct_;
            
            //! The object used to compute the search direction for every loop iteration.
            /*!
             *  By default, it is set equal to an object of type dcp::GradientSearchDirection defaul-constructed,
             *  but it can be changed using the function \c setSearchDirection
             */
             SearchDirectionComputer searchDirectionComputer_;

            // ---------------------------------------------------------------------------------------------//

        private:
            
    };
}

#endif
