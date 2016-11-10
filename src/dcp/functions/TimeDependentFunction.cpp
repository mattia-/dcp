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

#include <dolfin/math/basic.h>
#include <dcp/functions/TimeDependentFunction.h>
#include <dcp/expressions/TimeDependentExpression.h>

namespace dcp
{
    TimeDependentFunction::TimeDependentFunction () :
        std::vector<TimeDependentFunction::value_type> (),
        parameters ("time_dependent_solution_parameters")
    {}



    TimeDependentFunction::TimeDependentFunction (const int &n,
                                                  const double& time,
                                                  const std::shared_ptr<const dolfin::FunctionSpace> functionSpace) :
        std::vector<TimeDependentFunction::value_type> (n, std::make_pair (time, dolfin::Function (functionSpace))),
        parameters ("time_dependent_solution_parameters")
    {}



    TimeDependentFunction::TimeDependentFunction (const int &n,
                                                  const TimeDependentFunction::value_type& value) :
        std::vector<TimeDependentFunction::value_type> (n, value),
        parameters ("time_dependent_solution_parameters")
    {}


    TimeDependentFunction& TimeDependentFunction::operator= (dolfin::Expression& expression)
    {
        // try to convert pointer to expressin to dcp::TimeDependentExpression
        bool isTimeDependentExpression = false;
        auto pointerToExpression = dynamic_cast<dcp::TimeDependentExpression*> (&expression);
        if (pointerToExpression != nullptr)
        {
            isTimeDependentExpression = true;
        }

        for (auto& element : (*this))
        {
            if (isTimeDependentExpression)
            {
                pointerToExpression->time ()->setTo (element.first);
            }

            element.second = expression;
        }

        return *this;
    }




    // ****************************************** //
    // ********** FRIEND FUNCTIONS ************** //
    // ****************************************** //
    TimeDependentFunction operator+ (const TimeDependentFunction& left,
                                     const TimeDependentFunction& right)
    {
        if (left.size () != right.size ())
        {
            dolfin::dolfin_error ("dcp: TimeDependentFunction.cpp",
                                  "perform sum between objects of type dcp::TimeDependentFunction in operator+",
                                  "Size mismatch");
        }

        // the two function spaces are supposed to be equal; otherwise, dolfin will issue an error when trying to sum
        // the functions anyway, so no need to check for it now
        TimeDependentFunction result;
        result.reserve (left.size ());

        // loop through time steps and sum functions
        dolfin::begin (dolfin::DBG, "Summing time dependent functions...");
        for (auto i = 0; i < left.size (); ++i)
        {
            if (dolfin::near (left[i].first, right[i].first) == false)
            {
                dolfin::dolfin_error ("dcp: TimeDependentFunction.cpp",
                                      "perform sum between objects of type dcp::TimeDependentFunction in operator+",
                                      "Mismatch in time value at position %d; the two values are %f and %f",
                                      i,
                                      left[i].first,
                                      right[i].first);
            }

            // create function to hold the sum and store it into result
            dolfin::Function sum (left[i].second.function_space ());
            sum = left[i].second + right[i].second;
            result.emplace_back (std::make_pair (left[i].first, sum));
        }
        dolfin::end (); // Summing time dependent functions

        return result;
    }



    TimeDependentFunction operator* (const TimeDependentFunction& function,
                                     const double& scalar)
    {
        TimeDependentFunction result;
        result.reserve (function.size ());

        // loop through time steps and multiply functions by scalar
        dolfin::begin (dolfin::DBG, "Multipling time dependent function by scalar...");
        for (auto i = 0; i < function.size (); ++i)
        {
            // create function to hold the product and store it into result
            dolfin::Function product (function[i].second.function_space ());
            product = function[i].second * scalar;
            result.emplace_back (std::make_pair (function[i].first, product));
        }
        dolfin::end (); // Summing time dependent functions

        return result;
    }



    TimeDependentFunction operator* (const double& scalar,
                                     const TimeDependentFunction& function)
    {
        return function * scalar;
    }




    TimeDependentFunction operator* (const TimeDependentFunction& function,
                                     const int& scalar)
    {
        return function * double(scalar);
    }



    TimeDependentFunction operator* (const int& scalar,
                                     const TimeDependentFunction& function)
    {
        return function * double(scalar);
    }




    TimeDependentFunction operator- (const TimeDependentFunction& function)
    {
        return (-1) * function;
    }
}
