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
    TimeDependentExpression::TimeDependentExpression (const double& t) : 
        VariableExpression (),
        t_ (t)
    {
        dolfin::log (dolfin::DBG, "TimeDependentExpression object created");
    }
    


    TimeDependentExpression::TimeDependentExpression (std::size_t dim, const double& t) : 
        VariableExpression (dim),
        t_ (t)
    { 
        dolfin::log (dolfin::DBG, "TimeDependentExpression object created");
    }
    



    TimeDependentExpression::TimeDependentExpression (std::size_t dim0, std::size_t dim1, const double& t) : 
        VariableExpression (dim0, dim1),
        t_ (t)
    { 
        dolfin::log (dolfin::DBG, "TimeDependentExpression object created");
    }
    


    TimeDependentExpression::TimeDependentExpression (std::vector<std::size_t> value_shape, const double& t) : 
        VariableExpression (value_shape),
        t_ (t)
    { 
        dolfin::log (dolfin::DBG, "TimeDependentExpression object created");
    }

    

    TimeDependentExpression::
    TimeDependentExpression (const std::map <std::string, std::shared_ptr <const dolfin::GenericFunction>>& variables, 
                             const double& t) : 
        VariableExpression (),
        t_ (t)
    {
        dolfin::log (dolfin::DBG, "TimeDependentExpression object created");
    }
        
           

    TimeDependentExpression::
    TimeDependentExpression (std::size_t dim,
                             const std::map <std::string, std::shared_ptr <const dolfin::GenericFunction>>& variables, 
                             const double& t) :
        VariableExpression (dim),
        t_ (t)
    {
        dolfin::log (dolfin::DBG, "TimeDependentExpression object created");
    }

    

    TimeDependentExpression::
    TimeDependentExpression (std::size_t dim0, 
                             std::size_t dim1,
                             const std::map <std::string, std::shared_ptr <const dolfin::GenericFunction>>& variables,
                             const double& t) : 
        VariableExpression (dim0, dim1, variables),
        t_ (t)
    {
        dolfin::log (dolfin::DBG, "TimeDependentExpression object created");
    }

    

    TimeDependentExpression::
    TimeDependentExpression (std::vector<std::size_t> value_shape,
                             const std::map <std::string, std::shared_ptr <const dolfin::GenericFunction>>& variables, 
                             const double& t) : 
        VariableExpression (value_shape, variables),
        t_ (t)
    {
        dolfin::log (dolfin::DBG, "TimeDependentExpression object created");
    }
        
    

    /******************* GETTERS *******************/
    double& TimeDependentExpression::t ()
    {
        return t_;
    }

    

    const double& TimeDependentExpression::t () const
    {
        return t_;
    }

            
            
    /******************* METHODS *******************/
    void TimeDependentExpression::eval (dolfin::Array<double>& values, 
                                        const dolfin::Array<double>& x,
                                        const ufc::cell& cell) const
    {
        // Redirect to simple eval
        eval (values, x); 
    }



    void TimeDependentExpression::eval (dolfin::Array<double>& values, const dolfin::Array<double>& x) const
    {
        dolfin::dolfin_error("TimeDependentExpression.cpp",
                             "evaluate expression",
                             "Missing eval() function (must be overloaded)");
    }



    void TimeDependentExpression::evaluateVariable (const std::string& variableName, 
                                                    dolfin::Array<double>& values, 
                                                    const dolfin::Array<double>& x) const
    {
        auto variable = variables_.find (variableName);

        if (variable == variables_.end ())
        {
            dolfin::dolfin_error ("dcp: TimeDependentExpression.cpp",
                                  "evaluateVariable", 
                                  "Cannot find variable \"%s\" in TimeDependentExpression map", 
                                  variableName.c_str ());
        }

        (variable -> second) -> eval (values, x);
    }
}
