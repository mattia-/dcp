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
#include <dolfin/plot/plot.h>
#include <dcp/functions/TimeDependentFunction.h>
#include <dcp/expressions/TimeDependentExpression.h>

namespace dcp
{
    TimeDependentFunction::TimeDependentFunction () :
        std::vector<TimeDependentFunction::value_type> (),
        parameters ("time_dependent_solution_parameters"),
        plotter_ ()
    {}



    TimeDependentFunction::TimeDependentFunction (const int &n,
                                                  const double& t0,
                                                  const double& dt,
                                                  const std::shared_ptr<const dolfin::FunctionSpace> functionSpace) :
        std::vector<TimeDependentFunction::value_type> (),
        parameters ("time_dependent_solution_parameters"),
        plotter_ ()
    {
        this->reserve (n);

        for (auto i = 0; i < n; ++i)
        {
            this->emplace_back (std::make_pair (t0 + i * dt, dolfin::Function (functionSpace)));
        }
    }



    TimeDependentFunction::TimeDependentFunction (const double& t0,
                                                  const double& dt,
                                                  const double& tf,
                                                  const std::shared_ptr<const dolfin::FunctionSpace> functionSpace) :
        std::vector<TimeDependentFunction::value_type> (),
        parameters ("time_dependent_solution_parameters"),
        plotter_ ()
    {
        double currentTime = t0;
        std::size_t counter = 0;
        while (currentTime <= tf)
        {
            this->emplace_back (std::make_pair (currentTime, dolfin::Function (functionSpace)));

            counter++;
            currentTime = t0 + counter * dt;
        }
    }



    TimeDependentFunction::TimeDependentFunction (const int &n,
                                                  const TimeDependentFunction::value_type& value) :
        std::vector<TimeDependentFunction::value_type> (n, value),
        parameters ("time_dependent_solution_parameters")
    {}



    void TimeDependentFunction::plot (std::string title, const bool& pause)
    {
        if (title.empty ())
        {
            title = "Time = ";
        }
        else
        {
            title += ", time = ";
        }

        dolfin::begin (dolfin::DBG, "Plotting time-dependent function...");
        std::string titleWithTime;
        for (const auto& element : (*this))
        {
            titleWithTime = title + std::to_string (element.first);
            if (plotter_ == nullptr)
            {
                dolfin::log (dolfin::DBG, "Plotting in new dolfin::VTKPlotter object...");
                plotter_ = dolfin::plot (dolfin::reference_to_no_delete_pointer (element.second), titleWithTime);
            }
            else if (plotter_ -> is_compatible (dolfin::reference_to_no_delete_pointer (element.second)) == false)
            {
                dolfin::log (dolfin::DBG, "Existing plotter is not compatible with object to be plotted.");
                dolfin::log (dolfin::DBG, "Creating new dolfin::VTKPlotter object...");
                plotter_ = dolfin::plot (dolfin::reference_to_no_delete_pointer (element.second), titleWithTime);
            }
            else 
            {
                plotter_ -> parameters ["title"] = titleWithTime;
                plotter_ -> plot (dolfin::reference_to_no_delete_pointer (element.second));
            }
            
            if (pause)
            {
                dolfin::interactive ();
            }
        }
        dolfin::end (); // Plotting time-dependent function...
    }



    TimeDependentFunction& TimeDependentFunction::operator= (dolfin::Expression& expression)
    {
        // try to convert pointer to expressin to dcp::TimeDependentExpression
        bool isTimeDependentExpression = false;
        auto pointerToExpression = dynamic_cast<dcp::TimeDependentExpression*> (&expression);
        if (pointerToExpression != nullptr)
        {
            isTimeDependentExpression = true;
        }

        double initialTime = 0;
        if (isTimeDependentExpression)
        {
            initialTime = pointerToExpression->time ()->value ();
        }

        dolfin::begin (dolfin::DBG, "Setting time-dependent function values from expression...");
        for (auto& element : (*this))
        {
            dolfin::begin (dolfin::DBG, "Setting function values for time %f...", element.first);
            if (isTimeDependentExpression)
            {
                pointerToExpression->setTime (element.first);
            }

            // use expression here (not pointerToExpression), so that every kind of expression can be assigned to
            // elemen.second
            element.second = expression;
            dolfin::end (); // Setting function values for time %f...
        }
        dolfin::end (); // Setting time-dependent function values from expression

        if (isTimeDependentExpression)
        {
            pointerToExpression->setTime (initialTime);
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
        for (std::size_t i = 0; i < left.size (); ++i)
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



    TimeDependentFunction operator- (const TimeDependentFunction& left,
                                     const TimeDependentFunction& right)
    {
        return left + (-right);
    }



    TimeDependentFunction operator* (const TimeDependentFunction& function,
                                     const double& scalar)
    {
        TimeDependentFunction result;
        result.reserve (function.size ());

        // loop through time steps and multiply functions by scalar
        dolfin::begin (dolfin::DBG, "Multipling time dependent function by scalar...");
        for (std::size_t i = 0; i < function.size (); ++i)
        {
            // create function to hold the product and store it into result
            dolfin::Function product (function[i].second.function_space ());
            product = function[i].second * scalar;
            result.emplace_back (std::make_pair (function[i].first, product));
        }
        dolfin::end (); // Multipling time dependent functions

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



    TimeDependentFunction operator- (const dcp::TimeDependentFunction& function)
    {
        return function * (-1);
    }



    void operator<< (dolfin::File& file, const dcp::TimeDependentFunction& function)
    {
        try
        {
            for (const auto& element : function)
            {
                file << std::make_pair (&element.second, element.first);
            }
        }
        catch (std::runtime_error& e)
        {
            for (const auto& element : function)
            {
                file << element.second;
            }
        }
    }
}
