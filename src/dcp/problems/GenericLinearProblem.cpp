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

#include "GenericLinearProblem.h"

namespace dcp
{
    /************************* CONSTRUCTORS ********************/
    GenericLinearProblem::GenericLinearProblem (const std::shared_ptr<const dolfin::FunctionSpace> functionSpace) :
        GenericProblem (functionSpace)
    {
        dolfin::begin (dolfin::DBG, "Building GenericLinearProblem...");

        dolfin::log (dolfin::DBG, "Setting up parameters...");
        parameters.add ("problem_type", "linear");

        dolfin::end ();

        dolfin::log (dolfin::DBG, "GenericLinearProblem object created");
    }
}
