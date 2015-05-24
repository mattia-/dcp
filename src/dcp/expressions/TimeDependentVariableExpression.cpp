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

#include <dcp/expressions/TimeDependentVariableExpression.h>
#include <dcp/expressions/DefaultEvaluator.h>
#include <utility>
#include <memory>

namespace dcp
{
    /******************* CONSTRUCTORS *******************/
    TimeDependentVariableExpression::TimeDependentVariableExpression 
        (const dcp::TimeDependentVariableExpression::Evaluator& evaluator,
         std::shared_ptr<dcp::Time> time) 
        : 
        dcp::TimeDependentExpression (dcp::DefaultEvaluator (), time),
        evaluator_ (evaluator)
    {
        dolfin::log (dolfin::DBG, "TimeDependentVariableExpression object created");
    }



    TimeDependentVariableExpression::TimeDependentVariableExpression 
        (std::size_t dim, 
         const dcp::TimeDependentVariableExpression::Evaluator& evaluator,
         std::shared_ptr<dcp::Time> time) 
        : 
        dcp::TimeDependentExpression (dim, dcp::DefaultEvaluator (), time),
        evaluator_ (evaluator)
    { 
        dolfin::log (dolfin::DBG, "TimeDependentVariableExpression object created");
    }




    TimeDependentVariableExpression::TimeDependentVariableExpression 
        (std::size_t dim0, 
         std::size_t dim1,
         const dcp::TimeDependentVariableExpression::Evaluator& evaluator,
         std::shared_ptr<dcp::Time> time) 
        : 
        dcp::TimeDependentExpression (dim0, dim1, dcp::DefaultEvaluator (), time),
        evaluator_ (evaluator)
    { 
        dolfin::log (dolfin::DBG, "TimeDependentVariableExpression object created");
    }



    TimeDependentVariableExpression::TimeDependentVariableExpression 
        (std::vector<std::size_t> value_shape,
         const dcp::TimeDependentVariableExpression::Evaluator& evaluator,
         std::shared_ptr<dcp::Time> time) 
        : 
        dcp::TimeDependentExpression (value_shape, dcp::DefaultEvaluator (), time),
        evaluator_ (evaluator)
    { 
        dolfin::log (dolfin::DBG, "TimeDependentVariableExpression object created");
    }



    TimeDependentVariableExpression::TimeDependentVariableExpression 
        (const std::map <std::string, std::shared_ptr <const dolfin::GenericFunction>>& variables,
         const dcp::TimeDependentVariableExpression::Evaluator& evaluator,
         std::shared_ptr<dcp::Time> time) 
        : 
        dcp::TimeDependentExpression (variables, dcp::DefaultEvaluator (), time),
        evaluator_ (evaluator)
    {
        dolfin::log (dolfin::DBG, "TimeDependentVariableExpression object created");
    }



    TimeDependentVariableExpression::TimeDependentVariableExpression 
        (std::size_t dim,
         const std::map <std::string, std::shared_ptr <const dolfin::GenericFunction>>& variables,
         const dcp::TimeDependentVariableExpression::Evaluator& evaluator,
         std::shared_ptr<dcp::Time> time) 
        : 
        dcp::TimeDependentExpression (dim, variables, dcp::DefaultEvaluator (), time),
        evaluator_ (evaluator)
    {
        dolfin::log (dolfin::DBG, "TimeDependentVariableExpression object created");
    }



    TimeDependentVariableExpression::TimeDependentVariableExpression 
        (std::size_t dim0, 
         std::size_t dim1,
         const std::map <std::string, std::shared_ptr <const dolfin::GenericFunction>>& variables,
         const dcp::TimeDependentVariableExpression::Evaluator& evaluator,
         std::shared_ptr<dcp::Time> time)
        :
        dcp::TimeDependentExpression (dim0, dim1, variables, dcp::DefaultEvaluator (), time),
        evaluator_ (evaluator)
    {
        dolfin::log (dolfin::DBG, "TimeDependentVariableExpression object created");
    }



    TimeDependentVariableExpression::TimeDependentVariableExpression 
        (std::vector<std::size_t> value_shape,
         const std::map <std::string, std::shared_ptr <const dolfin::GenericFunction>>& variables,
         const dcp::TimeDependentVariableExpression::Evaluator& evaluator,
         std::shared_ptr<dcp::Time> time)
        :
        dcp::TimeDependentExpression (value_shape, variables, dcp::DefaultEvaluator (), time),
        evaluator_ (evaluator)
    {
        dolfin::log (dolfin::DBG, "TimeDependentVariableExpression object created");
    }



    /******************* METHODS *******************/
    void TimeDependentVariableExpression::eval (dolfin::Array<double>& values, 
                                                const dolfin::Array<double>& x, 
                                                const ufc::cell& cell) const
    {
        // redirect to simple eval
        this -> eval (values, x);
    }
    
    

    void TimeDependentVariableExpression::eval (dolfin::Array<double>& values, const dolfin::Array<double>& x) const
    {
        evaluator_ (values, x, time_ -> value (), variables_);
    }
    


    dcp::TimeDependentVariableExpression* TimeDependentVariableExpression::clone () const
    {
        return new dcp::TimeDependentVariableExpression (*this);
    }
}
