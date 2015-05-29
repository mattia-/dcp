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

#include <dcp/optimizers/AbstractDescentMethod.h>
#include <dolfin/log/dolfin_log.h>

namespace dcp
{
    AbstractDescentMethod::AbstractDescentMethod () : 
        parameters ("optimizer_parameters"),
        dotProduct_ (),
        searchDirectionComputer_ (dcp::GradientSearchDirection ())
    {
        dolfin::log (dolfin::DBG, "AbstractDescentMethod object created");
    }
    


    void AbstractDescentMethod::setDotProduct (const dolfin::Form& dotProductForm)
    {
        dotProduct_.setDotProductComputer (dotProductForm);
    }
    


    void AbstractDescentMethod::setSearchDirection (const SearchDirectionComputer& searchDirectionComputer)
    {
        searchDirectionComputer_ = searchDirectionComputer;
    }
}
