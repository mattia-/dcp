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

#ifndef SRC_EXPRESSIONS_TIMEDEPENDENTEXPRESSION_H_INCLUDE_GUARD
#define SRC_EXPRESSIONS_TIMEDEPENDENTEXPRESSION_H_INCLUDE_GUARD

#include <dolfin/function/GenericFunction.h>
#include <dolfin/function/Expression.h>
#include <dolfin/function/Function.h>
#include <dolfin/common/Array.h>
#include <map>
#include <vector>
#include <string>
#include <dcp/expressions/VariableExpression.h>

namespace dcp
{
    /*! class TimeDependentExpression TimeDependentExpression.h
     *  \brief Class for expressions depending on time. It derives from \c dcp::VariableExpression so it inherits
     *  all of its functionalities and can depend on other expression too.
     *  
     *  This class represents an expression which can depend on a time parameter as well as the value of other 
     *  expressions and functions. In general, a concrete class will need to be
     *  defined as deriving from this one and the method \c eval() will have to be overridden with a user-defined
     *  expression, much like what happens with \c dolfin::Expression itself. 
     *  The time parameter is stored in the protected member \c t_ and can be accessed through the methods
     *  \c t() and \c setTime(). Remember that the constructor does not accept the time as input argument because that
     *  would force the user to implement a constructor of his own when he derives from this class to implement the
     *  \c eval() method, so \c t_ will be default-constructed. Before using the class (and all of those that derive 
     *  from it), the method \c setTime() must be called if \c t_ should be different from 0.
     */
    class TimeDependentExpression : public dcp::VariableExpression
    {
        // ---------------------------------------------------------------------------------------------//  
        public:
            /******************* CONSTRUCTORS *******************/
            //! Default constructor. Create scalar expression
            TimeDependentExpression ();
            
            //! Create vector-valued expression with given dimension. This will call the appropriate 
            //! \c dolfin::Expression constructor
            /*
             *  Input arguments:
             *  \param dim dimension of the vector-valued expression
             */         
            explicit TimeDependentExpression (std::size_t dim);

            //! Create matrix-valued expression with given dimensions. This will call the appropriate 
            //! \c dolfin::Expression constructor
            /*!
             *  Input arguments:
             *  \param dim0 dimension (rows)
             *  \param dim1 dimension (columns)
             */          
            TimeDependentExpression (std::size_t dim0, std::size_t dim1);

            //! Create tensor-valued expression with given shape. This will call the appropriate \c dolfin::Expression
            //! constructor
            /*!
             *  Input arguments:
             *  \param value_shape shape of expression
             */          
            explicit TimeDependentExpression (std::vector<std::size_t> value_shape);

            //! Constructor from \c std::map
            /*!
             *  Uses \c map passed as input to create the protected member \c variables_
             *  Input arguments:
             *  \param variables map used to initialize the protected member \c variables_
             */
            TimeDependentExpression 
                (const std::map <std::string, std::shared_ptr <const dolfin::GenericFunction>>& variables);
            
            //! Create vector-valued expression with given dimension and given map. This will call the appropriate 
            //! \c dolfin::Expression constructor and set the protected member \c variables_ using the input \c map
            /*
             *  Input arguments:
             *  \param dim dimension of the vector-valued expression
             *  \param variables map used to initialize the protected member \c variables_
             */         
            explicit TimeDependentExpression 
                (std::size_t dim,
                 const std::map <std::string, std::shared_ptr <const dolfin::GenericFunction>>& variables);

            //! Create matrix-valued expression with given dimension and given map. This will call the appropriate 
            //! \c dolfin::Expression constructor and set the protected member \c variables_ using the input \c map
            /*
             *  Input arguments:
             *  \param dim0 dimension (rows)
             *  \param dim1 dimension (columns)
             *  \param variables map used to initialize the protected member \c variables_
             */         
            TimeDependentExpression 
                (std::size_t dim0, 
                 std::size_t dim1,
                 const std::map <std::string, std::shared_ptr <const dolfin::GenericFunction>>& variables);

            //! Create tensor-valued expression with given dimension and given map. This will call the appropriate 
            //! \c dolfin::Expression constructor and set the protected member \c variables_ using the input \c map
            /*
             *  Input arguments:
             *  \param value_shape shape of expression
             *  \param variables map used to initialize the protected member \c variables_
             */         
            explicit TimeDependentExpression
                (std::vector<std::size_t> value_shape,
                 const std::map <std::string, std::shared_ptr <const dolfin::GenericFunction>>& variables);

            //! Default copy constructor
            /*!
             *  Input arguments:
             *  \param expression object to be copied
             */
            TimeDependentExpression (const TimeDependentExpression& expression) = default;
            
            //! Default move constructor
            /*!
             *  Input arguments:
             *  \param expression object to be moved
             */
            TimeDependentExpression (TimeDependentExpression&& expression) = default;
            
            
            /******************* DESTRUCTOR *******************/
            //! Default destructor                
            /*! 
             *  Default destructor, since members of the class are trivially 
             *  destructible.
             */
            virtual ~TimeDependentExpression () {};


            /******************* GETTERS *******************/
            //! Get a const reference to the current time stored in the class
            const double& t () const;


            /******************* SETTERS *******************/
            //! Set the current time
            void setTime (const double& t);

            
            /******************* METHODS *******************/
            //! Evaluate at given point in given cell. Overrides method in \c dcp::VariableExpression
            /*!
             *  Input arguments are:
             *  \param values array that will contain the evaluated function at the given point
             *  \param x the coordinates of the point
             *  \param cell the cell which contains the given point
             */
            virtual void eval (dolfin::Array<double>& values, 
                               const dolfin::Array<double>& x, 
                               const ufc::cell& cell) const override;

            //! Evaluate at given point in given cell. Overrides method in \c dcp::VariableExpression
            /*!
             *  Input arguments are:
             *  \param values array that will contain the evaluated function at the given point
             *  \param x the coordinates of the point
             */
            virtual void eval (dolfin::Array<double>& values, const dolfin::Array<double>& x) const override;

            //! Evaluate variable identified by given name at given point. Overrides method in 
            //! \c dcp::VariableExpression
            /*!
             *  Input arguments are:
             *  \param variableName string to identify the variable we want to evaluate
             *  \param values array that will contain the evaluated function at the given point
             *  \param x the coordinates of the point at which evaluate the variable
             */
            virtual void evaluateVariable (const std::string& variableName, 
                                           dolfin::Array<double>& values, 
                                           const dolfin::Array<double>& x) const override;


        // ---------------------------------------------------------------------------------------------//  
        protected:
            //! The time at which the time dependent expression should be evaluated
            double t_;
            
    };
}

#endif
