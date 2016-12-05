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

#include <dcp/problems/AlgebraicProblem.h>

namespace dcp
{
    /******************* CONSTRUCTORS *******************/
    AlgebraicProblem::AlgebraicProblem (const std::shared_ptr<const dolfin::FunctionSpace> functionSpace,
                                        const std::shared_ptr<dcp::GenericExpression> expression) :
        GenericProblem (functionSpace),
        expression_ (expression)
    { 
        dolfin::begin (dolfin::DBG, "Building AlgebraicProblem...");

        solution_.emplace_back (std::make_pair (-1, dolfin::Function (functionSpace_)));
        
        dolfin::log (dolfin::DBG, "Setting up parameters...");
        parameters.add ("problem_type", "algebraic");

        dolfin::end (); // "Building AlgebraicProblem"
        
        dolfin::log (dolfin::DBG, "AlgebraicProblem object created");
    }



    AlgebraicProblem::AlgebraicProblem (const std::shared_ptr<const dolfin::FunctionSpace> functionSpace,
                                        const dcp::GenericExpression& expression) :
        GenericProblem (functionSpace),
        expression_ (expression.clone ())
    { 
        dolfin::begin (dolfin::DBG, "Building AlgebraicProblem...");

        solution_.emplace_back (std::make_pair (-1, dolfin::Function (functionSpace_)));
        
        dolfin::log (dolfin::DBG, "Setting up parameters...");
        parameters.add ("problem_type", "algebraic");
        
        dolfin::end (); // "Building AlgebraicProblem"

        dolfin::log (dolfin::DBG, "AlgebraicProblem object created");
    }



    /******************* GETTERS *******************/
    const dcp::GenericExpression& AlgebraicProblem::expression () const
    {
        return *expression_;
    }

    

    /******************* SETTERS *******************/
    void AlgebraicProblem::setCoefficient (const std::string& coefficientType, 
                                           const std::shared_ptr<const dolfin::GenericFunction> coefficientValue,
                                           const std::string& coefficientName)
    {
        if (coefficientType == "expression")
        {
            dolfin::log (dolfin::DBG, "Setting expression coefficient \"%s\"...", coefficientName.c_str ());
            expression_->setCoefficient (coefficientName, coefficientValue);
        }
        else
        {
            dolfin::warning ("Cannot set coefficient \"%s\" in algebraic problem. Coefficient type \"%s\" unknown",
                             coefficientName.c_str (),
                             coefficientType.c_str ());
        }

    }



    void AlgebraicProblem::setCoefficient (const std::string& coefficientType,
                                           const std::shared_ptr<const dolfin::GenericFunction> coefficientValue,
                                           const std::size_t& coefficientNumber)
    {
        // Do nothing. Just overload base class function to make class not pure-virtual
    }



    void AlgebraicProblem::setIntegrationSubdomain (const std::string& formType,
                                                    std::shared_ptr<const dolfin::MeshFunction<std::size_t>> meshFunction,
                                                    const dcp::SubdomainType& subdomainType)
    {
        // Do nothing. Just overload base class function to make class not pure-virtual
    }



    bool AlgebraicProblem::addDirichletBC (const dolfin::GenericFunction& condition, 
                                           const dolfin::SubDomain& boundary,
                                           std::string bcName)
    {
        // Do nothing
        return true;
    }



    bool AlgebraicProblem::addDirichletBC (const dolfin::GenericFunction& condition, 
                                           const dolfin::SubDomain& boundary,
                                           const std::size_t& component,
                                           std::string bcName)
    {
        // Do nothing
        return true;
    }



    bool AlgebraicProblem::addDirichletBC (std::shared_ptr<const dolfin::GenericFunction> condition, 
                                           std::shared_ptr<const dolfin::SubDomain> boundary,
                                           std::string bcName)
    {
        // Do nothing
        return true;
    }

    

    bool AlgebraicProblem::addDirichletBC (std::shared_ptr<const dolfin::GenericFunction> condition, 
                                           std::shared_ptr<const dolfin::SubDomain> boundary,
                                           const std::size_t& component,
                                           std::string bcName)
    {
        // Do nothing
        return true;
    }



    bool AlgebraicProblem::addDirichletBC (const dolfin::DirichletBC& dirichletCondition, 
                                           std::string bcName)
    {
        // Do nothing
        return true;
    }


    bool AlgebraicProblem::addDirichletBC (dolfin::DirichletBC&& dirichletCondition, 
                                           std::string bcName)
    {
        // Do nothing
        return true;
    }



    bool AlgebraicProblem::removeDirichletBC (const std::string& bcName)
    {
        // Do nothing
        return true;
    }



    /******************* METHODS *******************/
    void AlgebraicProblem::solve (const std::string& solveType)
    {
        if (solveType != "default" && solveType != "stash")
        {
            dolfin::dolfin_error ("dcp: AlgebraicProblem.h", 
                                  "solve",
                                  "Unknown solve type \"%s\" requested",
                                  solveType.c_str ());
        }

        dolfin::log (dolfin::DBG, "Solve type: %s", solveType.c_str ());

        if (solveType == "default")
        {
            solution_.back ().second = *(expression_);

            // set stashedSolution_ to be equal to the last computed solution, so that when solution() is called
            // from a subiterations loop it gets the right one. In the case of "default" solveType, indeed, the
            // solution returned should be the same for all solution types
            stashedSolution_ = solution_.back ().second;;
        }
        else if (solveType == "stash")
        {
            stashedSolution_ = *(expression_);
        }
    }



    dcp::AlgebraicProblem* AlgebraicProblem::clone () const
    {
        dolfin::begin (dolfin::DBG, "Cloning object...");

        std::string cloneMethod = parameters["clone_method"];

        dolfin::log (dolfin::DBG, "Clone method: %s", cloneMethod.c_str ());
        dolfin::log (dolfin::DBG, "Creating new object of type AlgebraicProblem...");

        // create new object
        dcp::AlgebraicProblem* clonedProblem = nullptr;
        if (cloneMethod == "shallow_clone")
        {
            clonedProblem = new dcp::AlgebraicProblem (this->functionSpace_, this->expression_); 
        }
        else if (cloneMethod == "deep_clone")
        {
            std::shared_ptr<dolfin::FunctionSpace> functionSpaceCopy (new dolfin::FunctionSpace (*functionSpace_));

            clonedProblem = new dcp::AlgebraicProblem (functionSpaceCopy, *(this->expression_));
        }
        else
        {
            dolfin::dolfin_error ("dcp: AlgebraicProblem.h",
                                  "clone",
                                  "Cannot clone linear differential problem. Unknown clone method: \"%s\"",
                                  cloneMethod.c_str ());
        }

        // clear parameters set of newly created object so that it can be populated by the parameters of the object
        // being created.
        dolfin::log (dolfin::DBG, "Copying parameters to new object...");
        clonedProblem->parameters.clear ();
        clonedProblem->parameters = this->parameters;

        // copy solution
        dolfin::log (dolfin::DBG, "Copying solution...");
        clonedProblem->solution_ = this->solution_;

        dolfin::end ();

        return clonedProblem;
    }
}
