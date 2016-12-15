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

#include <dcp/time/Time_.h>
#include <dolfin/log/dolfin_log.h>

namespace dcp
{
    /******************* CONSTRUCTORS *******************/
    Time::Time (const double& time) : time_ (time)
    {
        dolfin::begin (dolfin::DBG, "Building Time...");

        dolfin::log (dolfin::DBG, "Set time value to %f", time_);

        dolfin::end ();

        dolfin::log (dolfin::DBG, "Time object created");
    }



    /******************* METHODS *******************/
    void Time::setTo (const double& time)
    {
        time_ = time;
    }



    void Time::add (const double& timeIncrement)
    {
        time_ += timeIncrement;
    }



    const double& Time::value () const
    {
        return time_;
    }
}
