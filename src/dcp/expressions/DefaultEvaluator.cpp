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

#include <dcp/expressions/DefaultEvaluator.h>

namespace dcp
{
    /******************* METHODS *******************/
    void DefaultEvaluator::operator() (dolfin::Array<double>& values, const dolfin::Array<double>& x)
    {
        dolfin::dolfin_error 
            ("dcp: Expression.cpp --> dcp::DefaultEvaluator",
             "evaluate expression",
             "Missing eval() function: you must either derive from dcp::Expression and override the eval() method or pass a functor to the constructor");
    }
    


    void DefaultEvaluator::operator() 
        (dolfin::Array<double>& values, 
         const dolfin::Array<double>& x, 
         const std::map <std::string, std::shared_ptr<const dolfin::GenericFunction> >& variables)
    {
        dolfin::dolfin_error 
            ("dcp: VariableExpression.cpp --> dcp::DefaultEvaluator",
             "evaluate expression",
             "Missing eval() function: you must either derive from dcp::VariableExpression and override the eval() method or pass a functor to the constructor");
    }
    


    void DefaultEvaluator::operator() 
        (dolfin::Array<double>& values, 
         const dolfin::Array<double>& x, 
         const double& t,
         const std::map <std::string, std::shared_ptr<const dolfin::GenericFunction> >& variables)
    {
        dolfin::dolfin_error 
            ("dcp: TimeDependentExpression.cpp --> dcp::DefaultEvaluator",
             "evaluate expression",
             "Missing eval() function: you must either derive from dcp::TimeDependentExpression and override the eval() method or pass a functor to the constructor");
    }
}

