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

#ifndef SRC_EXPRESSIONS_EXPRESSION_H_INCLUDE_GUARD
#define SRC_EXPRESSIONS_EXPRESSION_H_INCLUDE_GUARD

#include <dolfin/function/GenericFunction.h>
#include <dolfin/function/Expression.h>
#include <dolfin/function/Function.h>
#include <dolfin/common/Array.h>
#include <map>
#include <vector>
#include <string>
#include <functional>
#include <dcp/expressions/DefaultEvaluator.h>

namespace dcp
{
    /*! class Expression Expression.h
     *  \brief Class for expression definition.
     *  
     *  This class offers an alernative to the dolfin way to define expressions. 
     *  Instead of having the users derive a class from \c dolfin::Expression and override the \c eval() method,
     *  this class stores a function wrapper as protected member, initialized upon building through the 
     *  object passed to the constructor. This function wrapper is called when the method \c eval() is called.
     */
    class Expression : public dolfin::Expression
    {
        // ---------------------------------------------------------------------------------------------//  
        public:
            typedef std::function <void (dolfin::Array<double>&, const dolfin::Array<double>&)> Evaluator;
            
            /******************* CONSTRUCTORS *******************/
            //! Default constructor. Create scalar expression. 
            /*
             *  Input arguments
             *  \param evaluator the evaluator to be used when calling the \c eval() method. If no evaluator is passed,
             *  the default one will be used (which will just issue a \c dolfin_error : the behaviour in this case is
             *  the same as the normal <tt>dolfin::Expression</tt>s)
             */
            Expression (const Evaluator& evaluator = dcp::DefaultEvaluator ());
            
            //! Create vector-valued expression with given dimension. This will call the appropriate 
            //! \c dolfin::Expression constructor
            /*
             *  Input arguments:
             *  \param dim dimension of the vector-valued expression
             *  \param evaluator the evaluator to be used when calling the \c eval() method. If no evaluator is passed,
             *  the default one will be used (which will just issue a \c dolfin_error : the behaviour in this case is
             *  the same as the normal <tt>dolfin::Expression</tt>s)
             */         
            explicit Expression (std::size_t dim, const Evaluator& evaluator = dcp::DefaultEvaluator ());

            //! Create matrix-valued expression with given dimensions. This will call the appropriate 
            //! \c dolfin::Expression constructor
            /*!
             *  Input arguments:
             *  \param dim0 dimension (rows)
             *  \param dim1 dimension (columns)
             *  \param evaluator the evaluator to be used when calling the \c eval() method. If no evaluator is passed,
             *  the default one will be used (which will just issue a \c dolfin_error : the behaviour in this case is
             *  the same as the normal <tt>dolfin::Expression</tt>s)
             */          
            Expression (std::size_t dim0, std::size_t dim1, const Evaluator& evaluator = dcp::DefaultEvaluator ());

            //! Create tensor-valued expression with given shape. This will call the appropriate \c dolfin::Expression
            //! constructor
            /*!
             *  Input arguments:
             *  \param value_shape shape of expression
             *  \param evaluator the evaluator to be used when calling the \c eval() method. If no evaluator is passed,
             *  the default one will be used (which will just issue a \c dolfin_error : the behaviour in this case is
             *  the same as the normal <tt>dolfin::Expression</tt>s)
             */          
            explicit Expression (std::vector<std::size_t> value_shape, 
                                 const Evaluator& evaluator = dcp::DefaultEvaluator ());

            //! Default copy constructor
            /*!
             *  Input arguments:
             *  \param expression object to be copied
             */
            Expression (const Expression& expression) = default;
            
            //! Default move constructor
            /*!
             *  Input arguments:
             *  \param expression object to be moved
             */
            Expression (Expression&& expression) = default;
            
            
            /******************* DESTRUCTOR *******************/
            //! Default destructor                
            /*! 
             *  Default destructor, since members of the class are trivially 
             *  destructible.
             */
            virtual ~Expression () {};


            /******************* METHODS *******************/
            //! Evaluate at given point in given cell. Overrides method in \c dolfin::Expression
            /*!
             *  Input arguments are:
             *  \param values array that will contain the evaluated function at the given point
             *  \param x the coordinates of the point
             *  \param cell the cell which contains the given point
             */
            virtual void eval (dolfin::Array<double>& values, 
                               const dolfin::Array<double>& x, 
                               const ufc::cell& cell) const override;

            //! Evaluate at given point in given cell. Overrides method in \c dolfin::Expression
            /*!
             *  Input arguments are:
             *  \param values array that will contain the evaluated function at the given point
             *  \param x the coordinates of the point
             */
            virtual void eval (dolfin::Array<double>& values, const dolfin::Array<double>& x) const override;

        // ---------------------------------------------------------------------------------------------//  
        protected:

        // ---------------------------------------------------------------------------------------------//  
        private:
            //! The evaluator to use when the \c eval() method is called. Made private so that it cannot be used in
            //! derived classes, since it would make no sense. Derived classes should define their own evaluator
            Evaluator evaluator_;
    };
}

#endif
