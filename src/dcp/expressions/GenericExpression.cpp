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

#include <dcp/expressions/GenericExpression.h>
#include <utility>
#include <memory>

namespace dcp
{
    /******************* CONSTRUCTORS *******************/
    GenericExpression::GenericExpression () :
        dolfin::Expression (),
        variables_ ()
    {
        dolfin::log (dolfin::DBG, "GenericExpression object created");
    }



    GenericExpression::GenericExpression (std::size_t dim) :
        dolfin::Expression (dim),
        variables_ ()
    { 
        dolfin::log (dolfin::DBG, "GenericExpression object created");
    }




    GenericExpression::GenericExpression (std::size_t dim0, std::size_t dim1) :
        dolfin::Expression (dim0, dim1),
        variables_ ()
    { 
        dolfin::log (dolfin::DBG, "GenericExpression object created");
    }



    GenericExpression::GenericExpression (std::vector<std::size_t> value_shape) :
        dolfin::Expression (value_shape),
        variables_ ()
    { 
        dolfin::log (dolfin::DBG, "GenericExpression object created");
    }



    GenericExpression::GenericExpression 
        (const std::map <std::string, std::shared_ptr <const dolfin::GenericFunction>>& variables) 
        :
            dolfin::Expression (),
            variables_ (variables)
    {
        dolfin::log (dolfin::DBG, "GenericExpression object created");
    }



    GenericExpression:: GenericExpression 
        (std::size_t dim,
         const std::map <std::string, std::shared_ptr <const dolfin::GenericFunction>>& variables)
        :
            dolfin::Expression (dim), 
            variables_ (variables)
    {
        dolfin::log (dolfin::DBG, "GenericExpression object created");
    }



    GenericExpression:: GenericExpression 
        (std::size_t dim0, 
         std::size_t dim1,
         const std::map <std::string, std::shared_ptr <const dolfin::GenericFunction>>& variables)
        : 
            dolfin::Expression (dim0, dim1),
            variables_ (variables)
    {
        dolfin::log (dolfin::DBG, "GenericExpression object created");
    }



    GenericExpression:: GenericExpression 
        (std::vector<std::size_t> value_shape,
         const std::map <std::string, std::shared_ptr <const dolfin::GenericFunction>>& variables)
        : 
            dolfin::Expression (value_shape),
            variables_ (variables)
    {
        dolfin::log (dolfin::DBG, "GenericExpression object created");
    }



    /******************* SETTERS *******************/
    void GenericExpression::setCoefficient (const std::string& variableName, 
                                            const std::shared_ptr <const dolfin::GenericFunction> value)
    {
        dolfin::log (dolfin::DBG, "Inserting variable in map...");

        variables_ [variableName] = value;
    }



    /******************* GETTERS *******************/
    const dolfin::Function& GenericExpression::function (const std::string& variableName) const
    {
        auto variable = variables_.find (variableName);
        return *(std::dynamic_pointer_cast<const dolfin::Function> (variable -> second));
    }



    const dolfin::Expression& GenericExpression::expression (const std::string& variableName) const
    {
        auto variable = variables_.find (variableName);
        return *(std::dynamic_pointer_cast<const dolfin::Expression> (variable -> second));
    }



    std::vector<std::string> GenericExpression::variablesNames () const
    {
        std::vector<std::string> names (variables_.size ());
        
        std::size_t i = 0;
        for (auto& nameFunctionPair : variables_)
        {
            names [i] = nameFunctionPair.first;
        }
        
        return names;
    }
    


    /******************* METHODS *******************/
    void GenericExpression::eval (dolfin::Array<double>& values, 
                                  const dolfin::Array<double>& x, 
                                  const ufc::cell& cell) const
    {
        // redirect to simple eval
        this -> eval (values, x);
    }
    


    void GenericExpression::evaluateVariable (const std::string& variableName,
                                              dolfin::Array<double>& values,
                                              const dolfin::Array<double>& x) const
    {
        auto variable = variables_.find (variableName);

        if (variable == variables_.end ())
        {
            dolfin::dolfin_error ("dcp: GenericExpression.cpp",
                                  "evaluateVariable", 
                                  "Cannot find variable \"%s\" in GenericExpression map", 
                                  variableName.c_str ());
        }

        (variable -> second) -> eval (values, x);
    }
}
