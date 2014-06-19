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

#ifndef SRC_OBJECTIVEFUNCTIONAL_VARIABLEEXPRESSION_HPP_INCLUDE_GUARD
#define SRC_OBJECTIVEFUNCTIONAL_VARIABLEEXPRESSION_HPP_INCLUDE_GUARD

#include <dolfin/function/GenericFunction.h>
#include <dolfin/function/Expression.h>
#include <dolfin/function/Function.h>
#include <dolfin/common/Array.h>
#include <map>
#include <vector>
#include <string>

namespace DCP
{
    /*! class VariableExpression VariableExpression.hpp
     *  \brief Class for expressions containing variable members, such as \c dolfin::Functions and other 
     *  \c dolfin::Expression s
     *  
     *  This class represents an expression (that is a function defined in terms of its cooridnates instead of 
     *  its value in the mesh nodes) which can also depend on the value of other expressions and functions. 
     *  It derives from \c dolfin::Expression and it extends its funcionalities by introducing a map that stores
     *  the variable functions to be used by the \c eval() function. In general, a concrete class will need to be
     *  defined as deriving from this one and the method \c eval() will have to be overridden with a user-defined
     *  expression, much like what happens with \c dolfin::Expression itself
     */
    class VariableExpression : public dolfin::Expression
    {
        // ---------------------------------------------------------------------------------------------//  
        public:
            /******************* CONSTRUCTORS *******************/
            //! Default constructor. Create scalar expression
            VariableExpression () = default;
            
            //! Create vector-valued expression with given dimension. This will call the appropriate \c dolfin::Expression 
            //! constructor
            /*
             *  Input arguments:
             *  \param dim dimension of the vector-valued expression
             */         
            explicit VariableExpression (std::size_t dim);

            //! Create matrix-valued expression with given dimensions. This will call the appropriate \c dolfin::Expression
            //! constructor
            /*!
             *  Input arguments:
             *  \param dim0 dimension (rows)
             *  \param dim1 dimension (columns)
             */          
            VariableExpression (std::size_t dim0, std::size_t dim1);

            //! Create tensor-valued expression with given shape. This will call the appropriate \c dolfin::Expression
            //! constructor
            /*!
             *  Input arguments:
             *  \param value_shape shape of expression
             */          
            explicit VariableExpression (std::vector<std::size_t> value_shape);

            //! Constructor from \c std::map
            /*!
             *  Uses \c map passed as input to create the protected member \c variables_
             *  Input arguments:
             *  \param variables map used to initialize the protected member \c variables_
             */
            VariableExpression (const std::map <std::string, std::shared_ptr <const dolfin::GenericFunction>>& variables);
            
            //! Create vector-valued expression with given dimension and given map. This will call the appropriate 
            //! \c dolfin::Expression constructor and set the protected member \c variables_ using the input \c map
            /*
             *  Input arguments:
             *  \param dim dimension of the vector-valued expression
             *  \param variables map used to initialize the protected member \c variables_
             */         
            explicit VariableExpression (std::size_t dim,
                                         const std::map <std::string, std::shared_ptr <const dolfin::GenericFunction>>& variables);

            //! Create matrix-valued expression with given dimension and given map. This will call the appropriate 
            //! \c dolfin::Expression constructor and set the protected member \c variables_ using the input \c map
            /*
             *  Input arguments:
             *  \param dim0 dimension (rows)
             *  \param dim1 dimension (columns)
             *  \param variables map used to initialize the protected member \c variables_
             */         
            VariableExpression (std::size_t dim0, 
                                std::size_t dim1,
                                const std::map <std::string, std::shared_ptr <const dolfin::GenericFunction>>& variables);

            //! Create tensor-valued expression with given dimension and given map. This will call the appropriate 
            //! \c dolfin::Expression constructor and set the protected member \c variables_ using the input \c map
            /*
             *  Input arguments:
             *  \param value_shape shape of expression
             *  \param variables map used to initialize the protected member \c variables_
             */         
            explicit VariableExpression (std::vector<std::size_t> value_shape,
                                         const std::map <std::string, std::shared_ptr <const dolfin::GenericFunction>>& variables);

            //! Default copy constructor
            /*!
             *  Input arguments:
             *  \param expression object to be copied
             */
            VariableExpression (const VariableExpression& expression) = default;
            
            //! Default move constructor
            /*!
             *  Input arguments:
             *  \param expression object to be moved
             */
            VariableExpression (VariableExpression&& expression) = default;
            
            
            /******************* DESTRUCTOR *******************/
            //! Default destructor                
            /*! 
             *  Default destructor, since members of the class are trivially 
             *  destructible.
             */
            virtual ~VariableExpression () {};


            /******************* SETTERS *******************/
            //! Sets value for the variable with given name
            /*! 
             *  Input arguments are:
             *  \param variableName string identifying the variable we want to set
             *  \param value the value of such variable, given as a shared pointer to a \c dolfin::GenericFunction
             *  
             *  The pair created by the two input arguments will be inserted in the protected member \c variables_
             */
            virtual void setCoefficient (const std::string& variableName, 
                                 const std::shared_ptr <const dolfin::GenericFunction> value);


            /******************* GETTERS *******************/
            //! Get variable with given name and return it as a <tt> const dolfin::Function& </tt>
            /*! 
             *  Input arguments are:
             *  \param variableName the name of the variable to be returned. The method will look for an object with key
             *  matching \c name in the map \c variables_
             *
             *  \return the function, if found, as a <tt> const dolfin::Function& </tt>
             */
            virtual const dolfin::Function& function (const std::string& variableName) const;

            //! Get variable with given name and return it as a <tt> const dolfin::Expression& </tt>
            /*! 
             *  Input arguments are:
             *  \param variableName the name of the variable to be returned. The method will look for an object with key
             *  matching \c name in the map \c variables_
             *
             *  \return the function, if found, as a <tt> const dolfin::Expression& </tt>
             */
            virtual const dolfin::Expression& expression (const std::string& variableName) const;


            /******************* METHODS *******************/
            //! Evaluate at given point in given cell. Overrides method in \c dolfin::Expression
            /*!
             *  Input arguments are:
             *  \param values array that will contain the evaluated function at the given point
             *  \param x the coordinates of the point
             *  \param cell the cell which contains the given point
             */
            virtual void eval (dolfin::Array<double>& values, const dolfin::Array<double>& x, const ufc::cell& cell) const;

            //! Evaluate at given point in given cell. Overrides method in \c dolfin::Expression
            /*!
             *  Input arguments are:
             *  \param values array that will contain the evaluated function at the given point
             *  \param x the coordinates of the point
             */
            virtual void eval (dolfin::Array<double>& values, const dolfin::Array<double>& x) const;

            //! Evaluate variable identified by given name at given point
            /*!
             *  Input arguments are:
             *  \param variableName string to identify the variable we want to evaluate
             *  \param values array that will contain the evaluated function at the given point
             *  \param x the coordinates of the point at which evaluate the variable
             */
            virtual void evaluateVariable (const std::string& variableName, 
                                           dolfin::Array<double>& values, 
                                           const dolfin::Array<double>& x) const;

        // ---------------------------------------------------------------------------------------------//  
        protected:
            //! The map that associates variables' names and values
            std::map <std::string, std::shared_ptr<const dolfin::GenericFunction> > variables_;
    };
}

#endif
