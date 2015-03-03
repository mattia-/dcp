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

#include <dcp/expressions/TimeDependentExpression.h>
#include <utility>
#include <memory>

namespace dcp
{
    /******************* CONSTRUCTORS *******************/
    TimeDependentExpression::TimeDependentExpression (const dcp::TimeDependentExpression::Evaluator& evaluator) : 
        dcp::GenericExpression (),
        t_ (),
        t (t_),
        evaluator_ (evaluator)
    {
        dolfin::log (dolfin::DBG, "TimeDependentExpression object created");
    }



    TimeDependentExpression::TimeDependentExpression (std::size_t dim, 
                                                      const dcp::TimeDependentExpression::Evaluator& evaluator) : 
        dcp::GenericExpression (dim),
        t_ (),
        t (t_),
        evaluator_ (evaluator)
    { 
        dolfin::log (dolfin::DBG, "TimeDependentExpression object created");
    }




    TimeDependentExpression::TimeDependentExpression (std::size_t dim0, 
                                                      std::size_t dim1,
                                                      const dcp::TimeDependentExpression::Evaluator& evaluator) : 
        dcp::GenericExpression (dim0, dim1),
        t_ (),
        t (t_),
        evaluator_ (evaluator)
    { 
        dolfin::log (dolfin::DBG, "TimeDependentExpression object created");
    }



    TimeDependentExpression::TimeDependentExpression (std::vector<std::size_t> value_shape,
                                                      const dcp::TimeDependentExpression::Evaluator& evaluator) : 
        dcp::GenericExpression (value_shape),
        t_ (),
        t (t_),
        evaluator_ (evaluator)
    { 
        dolfin::log (dolfin::DBG, "TimeDependentExpression object created");
    }



    TimeDependentExpression::TimeDependentExpression 
        (const std::map <std::string, std::shared_ptr <dolfin::GenericFunction>>& variables,
         const dcp::TimeDependentExpression::Evaluator& evaluator) 
        : 
        dcp::GenericExpression (variables),
        t_ (),
        t (t_),
        evaluator_ (evaluator)
    {
        dolfin::log (dolfin::DBG, "TimeDependentExpression object created");
    }



    TimeDependentExpression::TimeDependentExpression 
        (std::size_t dim,
         const std::map <std::string, std::shared_ptr <dolfin::GenericFunction>>& variables,
         const dcp::TimeDependentExpression::Evaluator& evaluator) 
        : 
        dcp::GenericExpression (dim, variables),
        t_ (),
        t (t_),
        evaluator_ (evaluator)
    {
        dolfin::log (dolfin::DBG, "TimeDependentExpression object created");
    }



    TimeDependentExpression::TimeDependentExpression 
        (std::size_t dim0, 
         std::size_t dim1,
         const std::map <std::string, std::shared_ptr <dolfin::GenericFunction>>& variables,
         const dcp::TimeDependentExpression::Evaluator& evaluator)
        :
        dcp::GenericExpression (dim0, dim1, variables),
        t_ (),
        t (t_)
    {
        dolfin::log (dolfin::DBG, "TimeDependentExpression object created");
    }



    TimeDependentExpression::TimeDependentExpression 
        (std::vector<std::size_t> value_shape,
         const std::map <std::string, std::shared_ptr <dolfin::GenericFunction>>& variables,
         const dcp::TimeDependentExpression::Evaluator& evaluator)
        :
        dcp::GenericExpression (value_shape, variables),
        t_ (),
        t (t_),
        evaluator_ (evaluator)
    {
        dolfin::log (dolfin::DBG, "TimeDependentExpression object created");
    }



    /******************* GETTERS *******************/
    const double& TimeDependentExpression::time () const
    {
        return t_;
    }



    /******************* SETTERS *******************/
    void TimeDependentExpression::setTime (const double& time, const std::string& variableName)
    {
        if (variableName.empty ())
        {
            t_ = time;
            return;
        }
        
        // if variableName is not empty, get the variable with given name
        auto nameVariablePair = variables_.find (variableName);
        
        // if element was found and it is an object of type dcp::TimeDependentExpression (or derived from it),
        // call setTime()
        if (nameVariablePair != variables_.end () 
            &&
            std::dynamic_pointer_cast<dcp::TimeDependentExpression> (nameVariablePair->second) != nullptr)
        {
            std::dynamic_pointer_cast<dcp::TimeDependentExpression> (nameVariablePair->second) -> setTime (t);
        }
    } 



    void TimeDependentExpression::setTimeAll (const double& time)
    {
        // set time for this object
        setTime (time);
        
        // set time for all the variables
        for (auto nameVariablePair : variables_)
        {
            setTime (time, nameVariablePair.first);
        }
    }
            


    /******************* METHODS *******************/
    void TimeDependentExpression::eval (dolfin::Array<double>& values, const dolfin::Array<double>& x) const
    {
        evaluator_ (values, x, t_);
    }



    dcp::TimeDependentExpression* TimeDependentExpression::clone () const
    {
        return new dcp::TimeDependentExpression (*this);
    }
}
