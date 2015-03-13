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

#include <dcp/expressions/VariableExpression.h>
#include <utility>
#include <memory>

namespace dcp
{
    /******************* CONSTRUCTORS *******************/
    VariableExpression::VariableExpression (const dcp::VariableExpression::Evaluator& evaluator) :
        dcp::Expression (),
        evaluator_ (evaluator)
    {
        dolfin::log (dolfin::DBG, "VariableExpression object created");
    }



    VariableExpression::VariableExpression (std::size_t dim, 
                                            const dcp::VariableExpression::Evaluator& evaluator) :
        dcp::Expression (dim),
        evaluator_ (evaluator)
    { 
        dolfin::log (dolfin::DBG, "VariableExpression object created");
    }




    VariableExpression::VariableExpression (std::size_t dim0, std::size_t dim1, 
                                            const dcp::VariableExpression::Evaluator& evaluator) :
        dcp::Expression (dim0, dim1),
        evaluator_ (evaluator)
    { 
        dolfin::log (dolfin::DBG, "VariableExpression object created");
    }



    VariableExpression::VariableExpression (std::vector<std::size_t> value_shape, 
                                            const dcp::VariableExpression::Evaluator& evaluator) :
        dcp::Expression (value_shape),
        evaluator_ (evaluator)
    { 
        dolfin::log (dolfin::DBG, "VariableExpression object created");
    }



    VariableExpression::VariableExpression 
        (const std::map <std::string, std::shared_ptr <dolfin::GenericFunction>>& variables,
         const dcp::VariableExpression::Evaluator& evaluator) 
        : 
            dcp::Expression (variables),
            evaluator_ (evaluator)
    {
        dolfin::log (dolfin::DBG, "VariableExpression object created");
    }



    VariableExpression:: VariableExpression 
        (std::size_t dim,
         const std::map <std::string, std::shared_ptr <dolfin::GenericFunction>>& variables,
         const dcp::VariableExpression::Evaluator& evaluator) 
        :
            dcp::Expression (dim, variables), 
            evaluator_ (evaluator)
    {
        dolfin::log (dolfin::DBG, "VariableExpression object created");
    }



    VariableExpression:: VariableExpression 
        (std::size_t dim0, 
         std::size_t dim1,
         const std::map <std::string, std::shared_ptr <dolfin::GenericFunction>>& variables,
         const dcp::VariableExpression::Evaluator& evaluator) 
        : 
            dcp::Expression (dim0, dim1, variables),
            evaluator_ (evaluator)
    {
        dolfin::log (dolfin::DBG, "VariableExpression object created");
    }



    VariableExpression:: VariableExpression 
        (std::vector<std::size_t> value_shape,
         const std::map <std::string, std::shared_ptr <dolfin::GenericFunction>>& variables,
         const dcp::VariableExpression::Evaluator& evaluator) 
        : 
            dcp::Expression (value_shape, variables),
            evaluator_ (evaluator)
    {
        dolfin::log (dolfin::DBG, "VariableExpression object created");
    }



    /******************* METHODS *******************/
    void VariableExpression::eval (dolfin::Array<double>& values, const dolfin::Array<double>& x) const
    {
        evaluator_ (values, x, variables_);
    }



    dcp::VariableExpression* VariableExpression::clone () const
    {
        return new dcp::VariableExpression (*this);
    }
}
