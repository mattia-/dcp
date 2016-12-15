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

#include <dcp/expressions/Expression.h>
#include <utility>
#include <memory>

namespace dcp
{
    /******************* CONSTRUCTORS *******************/
    Expression::Expression (const dcp::Expression::Evaluator& evaluator) :
        dcp::GenericExpression (),
        evaluator_ (evaluator)
    {
        dolfin::log (dolfin::DBG, "Expression object created");
    }



    Expression::Expression (std::size_t dim, const dcp::Expression::Evaluator& evaluator) :
        dcp::GenericExpression (dim),
        evaluator_ (evaluator)
    {
        dolfin::log (dolfin::DBG, "Expression object created");
    }




    Expression::Expression (std::size_t dim0, std::size_t dim1, const dcp::Expression::Evaluator& evaluator) :
        dcp::GenericExpression (dim0, dim1),
        evaluator_ (evaluator)
    {
        dolfin::log (dolfin::DBG, "Expression object created");
    }



    Expression::Expression (std::vector<std::size_t> value_shape, const dcp::Expression::Evaluator& evaluator) :
        dcp::GenericExpression (value_shape),
        evaluator_ (evaluator)
    {
        dolfin::log (dolfin::DBG, "Expression object created");
    }



    Expression::Expression (const std::map <std::string, std::shared_ptr <const dolfin::GenericFunction>>& variables,
                            const dcp::Expression::Evaluator& evaluator) :
        dcp::GenericExpression (variables),
        evaluator_ (evaluator)
    {
        dolfin::log (dolfin::DBG, "Expression object created");
    }



    Expression:: Expression (std::size_t dim,
                             const std::map <std::string, std::shared_ptr <const dolfin::GenericFunction>>& variables,
                             const dcp::Expression::Evaluator& evaluator) :
        dcp::GenericExpression (dim, variables),
        evaluator_ (evaluator)
    {
        dolfin::log (dolfin::DBG, "Expression object created");
    }



    Expression:: Expression (std::size_t dim0,
                             std::size_t dim1,
                             const std::map <std::string, std::shared_ptr <const dolfin::GenericFunction>>& variables,
                             const dcp::Expression::Evaluator& evaluator) :
        dcp::GenericExpression (dim0, dim1, variables),
        evaluator_ (evaluator)
    {
        dolfin::log (dolfin::DBG, "Expression object created");
    }



    Expression:: Expression (std::vector<std::size_t> value_shape,
                             const std::map <std::string, std::shared_ptr <const dolfin::GenericFunction>>& variables,
                             const dcp::Expression::Evaluator& evaluator) :
        dcp::GenericExpression (value_shape, variables),
        evaluator_ (evaluator)
    {
        dolfin::log (dolfin::DBG, "Expression object created");
    }



    /******************* METHODS *******************/
    void Expression::eval (dolfin::Array<double>& values, const dolfin::Array<double>& x) const
    {
        evaluator_ (values, x);
    }



    dcp::Expression* Expression::clone () const
    {
        return new dcp::Expression (*this);
    }
}
