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
    TimeDependentExpression::TimeDependentExpression (const dcp::TimeDependentExpression::Evaluator& evaluator,
                                                      std::shared_ptr<dcp::Time> time) :
        dcp::GenericExpression (),
        time_ (time != nullptr ? time : std::shared_ptr<dcp::Time> (new dcp::Time (0))),
        evaluator_ (evaluator)
    {
        dolfin::log (dolfin::DBG, "TimeDependentExpression object created");
    }



    TimeDependentExpression::TimeDependentExpression (std::size_t dim,
                                                      const dcp::TimeDependentExpression::Evaluator& evaluator,
                                                      std::shared_ptr<dcp::Time> time) :
        dcp::GenericExpression (dim),
        time_ (time != nullptr ? time : std::shared_ptr<dcp::Time> (new dcp::Time (0))),
        evaluator_ (evaluator)
    {
        dolfin::log (dolfin::DBG, "TimeDependentExpression object created");
    }




    TimeDependentExpression::TimeDependentExpression (std::size_t dim0,
                                                      std::size_t dim1,
                                                      const dcp::TimeDependentExpression::Evaluator& evaluator,
                                                      std::shared_ptr<dcp::Time> time) :
        dcp::GenericExpression (dim0, dim1),
        time_ (time != nullptr ? time : std::shared_ptr<dcp::Time> (new dcp::Time (0))),
        evaluator_ (evaluator)
    {
        dolfin::log (dolfin::DBG, "TimeDependentExpression object created");
    }



    TimeDependentExpression::TimeDependentExpression (std::vector<std::size_t> value_shape,
                                                      const dcp::TimeDependentExpression::Evaluator& evaluator,
                                                      std::shared_ptr<dcp::Time> time) :
        dcp::GenericExpression (value_shape),
        time_ (time != nullptr ? time : std::shared_ptr<dcp::Time> (new dcp::Time (0))),
        evaluator_ (evaluator)
    {
        dolfin::log (dolfin::DBG, "TimeDependentExpression object created");
    }



    TimeDependentExpression::TimeDependentExpression
        (const std::map <std::string, std::shared_ptr <const dolfin::GenericFunction>>& variables,
         const dcp::TimeDependentExpression::Evaluator& evaluator,
         std::shared_ptr<dcp::Time> time)
        :
        dcp::GenericExpression (variables),
        time_ (time != nullptr ? time : std::shared_ptr<dcp::Time> (new dcp::Time (0))),
        evaluator_ (evaluator)
    {
        dolfin::log (dolfin::DBG, "TimeDependentExpression object created");
    }



    TimeDependentExpression::TimeDependentExpression
        (std::size_t dim,
         const std::map <std::string, std::shared_ptr <const dolfin::GenericFunction>>& variables,
         const dcp::TimeDependentExpression::Evaluator& evaluator,
         std::shared_ptr<dcp::Time> time)
        :
        dcp::GenericExpression (dim, variables),
        time_ (time != nullptr ? time : std::shared_ptr<dcp::Time> (new dcp::Time (0))),
        evaluator_ (evaluator)
    {
        dolfin::log (dolfin::DBG, "TimeDependentExpression object created");
    }



    TimeDependentExpression::TimeDependentExpression
        (std::size_t dim0,
         std::size_t dim1,
         const std::map <std::string, std::shared_ptr <const dolfin::GenericFunction>>& variables,
         const dcp::TimeDependentExpression::Evaluator& evaluator,
         std::shared_ptr<dcp::Time> time)
        :
        dcp::GenericExpression (dim0, dim1, variables),
        time_ (time != nullptr ? time : std::shared_ptr<dcp::Time> (new dcp::Time (0)))
    {
        dolfin::log (dolfin::DBG, "TimeDependentExpression object created");
    }



    TimeDependentExpression::TimeDependentExpression
        (std::vector<std::size_t> value_shape,
         const std::map <std::string, std::shared_ptr <const dolfin::GenericFunction>>& variables,
         const dcp::TimeDependentExpression::Evaluator& evaluator,
         std::shared_ptr<dcp::Time> time)
        :
        dcp::GenericExpression (value_shape, variables),
        time_ (time != nullptr ? time : std::shared_ptr<dcp::Time> (new dcp::Time (0))),
        evaluator_ (evaluator)
    {
        dolfin::log (dolfin::DBG, "TimeDependentExpression object created");
    }



    /******************* GETTERS *******************/
    std::shared_ptr<dcp::Time> TimeDependentExpression::time ()
    {
        return time_;
    }



    /******************* METHODS *******************/
    void TimeDependentExpression::setTime (const double& time)
    {
        time_ -> setTo (time);
    }



    void TimeDependentExpression::eval (dolfin::Array<double>& values, const dolfin::Array<double>& x) const
    {
        evaluator_ (values, x, time_ -> value ());
    }



    dcp::TimeDependentExpression* TimeDependentExpression::clone () const
    {
        return new dcp::TimeDependentExpression (*this);
    }
}
