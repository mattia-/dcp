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
#include <dcp/expressions/GenericExpression.h>

namespace dcp
{
    /*! class Expression Expression.h
     *  \brief Class for expression definition.
     *  
     *  This class offers an alernative to the dolfin way to define expressions. 
     *  It can be used in two ways:
     *  1) it can be derived from, exactly like one does for <tt>\c dolfin::Expression></tt>: concrete class will need 
     *  to be defined as deriving from this one and the method \c eval() will have to be overridden with a user-defined
     *  expression. 
     *  2) it can be constructed directly, passing a \c std::functional to the constructor. This functional will be 
     *  called when the \c eval() method is called.
     *  The former method is provided for ease of use, but the latter one should be preferred. In particular, note that
     *  if you choose to use the first method, you will probably want to override the \c clone() method as well in
     *  the derived class.
     */
    class Expression : public dcp::GenericExpression
    {
        // ---------------------------------------------------------------------------------------------//  
        public:
            typedef std::function <void (dolfin::Array<double>&, const dolfin::Array<double>&)> Evaluator;

            /******************* CONSTRUCTORS *******************/
            //! Default constructor. Create scalar expression. 
            /*
             *  Input arguments
             *  \param evaluator the evaluator to be used when calling the \c eval() method.
             */
            Expression (const dcp::Expression::Evaluator& evaluator);

            //! Create vector-valued expression with given dimension. This will call the appropriate 
            //! \c dolfin::Expression constructor
            /*
             *  Input arguments:
             *  \param dim dimension of the vector-valued expression
             *  \param evaluator the evaluator to be used when calling the \c eval() method.
             */         
            explicit Expression (std::size_t dim, const dcp::Expression::Evaluator& evaluator);

            //! Create matrix-valued expression with given dimensions. This will call the appropriate 
            //! \c dolfin::Expression constructor
            /*!
             *  Input arguments:
             *  \param dim0 dimension (rows)
             *  \param dim1 dimension (columns)
             *  \param evaluator the evaluator to be used when calling the \c eval() method.
             */          
            Expression (std::size_t dim0, std::size_t dim1, const dcp::Expression::Evaluator& evaluator);

            //! Create tensor-valued expression with given shape. This will call the appropriate \c dolfin::Expression
            //! constructor
            /*!
             *  Input arguments:
             *  \param value_shape shape of expression
             *  \param evaluator the evaluator to be used when calling the \c eval() method.
             */          
            explicit Expression (std::vector<std::size_t> value_shape, 
                                 const dcp::Expression::Evaluator& evaluator);

            //! Constructor from \c std::map
            /*!
             *  Uses \c map passed as input to create the protected member \c variables_
             *  Input arguments:
             *  \param variables map used to initialize the protected member \c variables_
             *  \param evaluator the evaluator to be used when calling the \c eval() method.
             */
            Expression (const std::map <std::string, std::shared_ptr <const dolfin::GenericFunction>>& variables,
                        const dcp::Expression::Evaluator& evaluator);

            //! Create vector-valued expression with given dimension and given map. This will call the appropriate 
            //! \c dolfin::Expression constructor and set the protected member \c variables_ using the input \c map
            /*
             *  Input arguments:
             *  \param dim dimension of the vector-valued expression
             *  \param variables map used to initialize the protected member \c variables_
             *  \param evaluator the evaluator to be used when calling the \c eval() method.
             */         
            explicit Expression (std::size_t dim,
                                 const std::map <std::string, std::shared_ptr <const dolfin::GenericFunction>>& variables,
                                 const dcp::Expression::Evaluator& evaluator);

            //! Create matrix-valued expression with given dimension and given map. This will call the appropriate 
            //! \c dolfin::Expression constructor and set the protected member \c variables_ using the input \c map
            /*
             *  Input arguments:
             *  \param dim0 dimension (rows)
             *  \param dim1 dimension (columns)
             *  \param variables map used to initialize the protected member \c variables_
             *  \param evaluator the evaluator to be used when calling the \c eval() method.
             */         
            Expression (std::size_t dim0, 
                        std::size_t dim1,
                        const std::map <std::string, std::shared_ptr <const dolfin::GenericFunction>>& variables,
                        const dcp::Expression::Evaluator& evaluator);


            //! Create tensor-valued expression with given dimension and given map. This will call the appropriate 
            //! \c dolfin::Expression constructor and set the protected member \c variables_ using the input \c map
            /*
             *  Input arguments:
             *  \param value_shape shape of expression
             *  \param variables map used to initialize the protected member \c variables_
             *  \param evaluator the evaluator to be used when calling the \c eval() method.
             */         
            explicit Expression (std::vector<std::size_t> value_shape,
                                 const std::map <std::string, std::shared_ptr <const dolfin::GenericFunction>>& variables,
                                 const dcp::Expression::Evaluator& evaluator);

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
            //! Evaluate at given point in given cell
            /*!
             *  Input arguments are:
             *  \param values array that will contain the evaluated function at the given point
             *  \param x the coordinates of the point
             */
            virtual void eval (dolfin::Array<double>& values, const dolfin::Array<double>& x) const override;

            //! Clone method
            /*!
             *  Performs a shallow clone
             *
             *  \return a pointer to the cloned object
             */
            virtual dcp::Expression* clone () const override;

            // ---------------------------------------------------------------------------------------------//  
        protected:

            // ---------------------------------------------------------------------------------------------//  
        private:
            //! The evaluator to use when the \c eval() method is called. Made private so that it cannot be used in
            //! derived classes, since it would make no sense. Derived classes should define their own evaluator
            dcp::Expression::Evaluator evaluator_;
    };
}

#endif
