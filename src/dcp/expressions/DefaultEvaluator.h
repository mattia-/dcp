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

#ifndef SRC_EXPRESSIONS_DEFAULTEVALUATOR_H_INCLUDE_GUARD
#define SRC_EXPRESSIONS_DEFAULTEVALUATOR_H_INCLUDE_GUARD

#include <dolfin/common/Array.h>
#include <dolfin/function/GenericFunction.h>
#include <string>
#include <map>
#include <memory>

namespace dcp
{
    /*! class DefaultEvaluator DefaultEvaluator.h
     *  \brief The default evaluator to be used in \c dcp::VariableExpression
     *  
     *  This class is just a functor whose call operator issues a \c dolfin_error, to mimic
     *  the behaviour of a \c dolfin::Expression when a \c dcp::VariableExpression is neither 
     *  derived from nor is built passing a functor to its constructor
     *  
     */
    class DefaultEvaluator
    {
        // ---------------------------------------------------------------------------------------------//  
        public:
            /******************* CONSTRUCTORS *******************/
            //! Default constructor. Create scalar expression
            DefaultEvaluator () = default;
            
            /******************* DESTRUCTOR *******************/
            //! Default destructor                
            /*! 
             *  Default destructor, since members of the class are trivially 
             *  destructible.
             */
            ~DefaultEvaluator () = default;

            /******************* METHODS *******************/
            //! Call operator for \c dcp::VariableExpression default evaluator
            void operator() (dolfin::Array<double>& values, 
                             const dolfin::Array<double>& x, 
                             const std::map <std::string, std::shared_ptr<const dolfin::GenericFunction> >& variables);
            
            
            //! Call operator for \c dcp::TimeDependentExpression default evaluator
            void operator() (dolfin::Array<double>& values, 
                             const dolfin::Array<double>& x, 
                             const double& t,
                             const std::map <std::string, std::shared_ptr<const dolfin::GenericFunction> >& variables);
            
        // ---------------------------------------------------------------------------------------------//  
        protected:
    };
}

#endif

