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

#include <splitting_methods/AbstractSplittingMethod.h>
#include <dolfin/log/dolfin_log.h>
#include <map>
#include <string>
#include <utility>

namespace dcp
{
    /************************* CONSTRUCTORS ********************/
    AbstractSplittingMethod::AbstractSplittingMethod (std::initializer_list<dolfin::FunctionSpace> functionSpaces) :
        parameters ("splitting_method_parameters"),
        functionSpaces_ (),
        differentialSystem_ (nullptr)
    { 
        dolfin::begin (dolfin::DBG, "Building AbstractSplittingMethod...");
        
        for (auto& functionSpace : functionSpaces)
        {
            functionSpaces_.push_back (dolfin::reference_to_no_delete_pointer (const_cast<dolfin::FunctionSpace&> (functionSpace)));
        }
            
        dolfin::end ();
        
        dolfin::log (dolfin::DBG, "AbstractSplittingMethod object created");
    }


    
    /************************* DESTRUCTOR ********************/
    AbstractSplittingMethod::~AbstractSplittingMethod ()
    {
        
    }

            

    /********************** GETTERS ***********************/
    const dcp::AbstractEquationSystem& AbstractSplittingMethod::system () const
    {
        return *differentialSystem_;
    }



    dcp::AbstractEquationSystem& AbstractSplittingMethod::system ()
    {
        return *differentialSystem_;
    }


    
    const dcp::AbstractProblem& AbstractSplittingMethod::operator[] (const std::string& name) const
    {
        return (*differentialSystem_) [name];
    }
    


    dcp::AbstractProblem& AbstractSplittingMethod::operator[] (const std::string& name)
    {
        return (*differentialSystem_) [name];
    }
    


    /********************** SETTERS ***********************/
    void AbstractSplittingMethod::setCoefficient (const std::string& problemName,
                                                  const std::string& coefficientType, 
                                                  const std::shared_ptr<const dolfin::GenericFunction> coefficientValue,
                                                  const std::string& coefficientName)
    {
        return (*differentialSystem_) [problemName].setCoefficient (coefficientType, coefficientValue, coefficientName);
    }
    


    void AbstractSplittingMethod::setCoefficient (const std::string& problemName,
                                                  const std::string& coefficientType,
                                                  const std::shared_ptr<const dolfin::GenericFunction> coefficientValue,
                                                  const std::size_t& coefficientNumber)
    {
        return (*differentialSystem_) [problemName].setCoefficient (coefficientType, coefficientValue, coefficientNumber);
    }
    


    void AbstractSplittingMethod::setIntegrationSubdomains 
        (const std::string& problemName,
         const std::string& formType,
         std::shared_ptr<const dolfin::MeshFunction<std::size_t>> meshFunction,
         const dcp::SubdomainType& subdomainType)
    {
        return (*differentialSystem_) [problemName].setIntegrationSubdomains (formType, meshFunction, subdomainType);
    }
    


    bool AbstractSplittingMethod::addDirichletBC (const std::string& problemName,
                                                  const dolfin::GenericFunction& condition, 
                                                  const dolfin::SubDomain& boundary,
                                                  std::string bcName)
    {
        return (*differentialSystem_) [problemName].addDirichletBC (condition, boundary, bcName); 
    }
    


    bool AbstractSplittingMethod::addDirichletBC (const std::string& problemName,
                                                  const dolfin::GenericFunction& condition, 
                                                  const dolfin::SubDomain& boundary, 
                                                  const std::size_t& component,
                                                  std::string bcName)
    {
        return (*differentialSystem_) [problemName].addDirichletBC (condition, boundary, component, bcName);
    }
    


    bool AbstractSplittingMethod::addDirichletBC (const std::string& problemName,
                                                  std::shared_ptr<const dolfin::GenericFunction> condition, 
                                                  std::shared_ptr<const dolfin::SubDomain> boundary,
                                                  std::string bcName)
    {
        return (*differentialSystem_) [problemName].addDirichletBC (condition, boundary, bcName);
    }

    

    bool AbstractSplittingMethod::addDirichletBC (const std::string& problemName,
                                                  std::shared_ptr<const dolfin::GenericFunction> condition, 
                                                  std::shared_ptr<const dolfin::SubDomain> boundary,
                                                  const std::size_t& component,
                                                  std::string bcName)
    {
        return (*differentialSystem_) [problemName].addDirichletBC (condition, boundary, component, bcName);
    }

    

    bool AbstractSplittingMethod::addDirichletBC (const std::string& problemName,
                                                  const dolfin::DirichletBC& dirichletCondition, 
                                                  std::string bcName)
    {
        return (*differentialSystem_) [problemName].addDirichletBC (dirichletCondition, bcName);
    }



    bool AbstractSplittingMethod::addDirichletBC (const std::string& problemName,
                                                  dolfin::DirichletBC&& dirichletCondition,
                                                  std::string bcName)
    {
        return (*differentialSystem_) [problemName].addDirichletBC (dirichletCondition, bcName);
    }



    bool AbstractSplittingMethod::removeDirichletBC (const std::string& problemName, const std::string& bcName)
    {
        return (*differentialSystem_) [problemName].removeDirichletBC (bcName);
    }
    
    

    /********************** METHODS ***********************/
    void AbstractSplittingMethod::apply (const std::string& type)
    {
        differentialSystem_ -> solve (true);
    } 
}
