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

#include <DifferentialProblem/TimeDependentDifferentialProblem.hpp>
#include <dolfin/log/dolfin_log.h>

namespace dcp
{
    /******************* CONSTRUCTORS *******************/

    TimeDependentDifferentialProblem::TimeDependentDifferentialProblem 
            (const std::shared_ptr<dolfin::Mesh> mesh, 
             const std::shared_ptr<dolfin::FunctionSpace> functionSpace,
             const double& dt,
             const double& endTime,
             const std::vector<std::string>& dtCoefficientTypes,
             const std::vector<std::string>& previousSolutionCoefficientTypes,
             const int& storeInterval,
             const std::shared_ptr<dcp::AbstractDifferentialProblem> timeSteppingProblem,
             const std::string& dtName,
             const std::string& previousSolutionName) :
        AbstractDifferentialProblem (mesh, functionSpace),
        timeSteppingProblem_ (timeSteppingProblem),
        solutions_ ()
    { 
        dolfin::begin (dolfin::DBG, "Building TimeDependentDifferentialProblem...");
        
        dolfin::log (dolfin::DBG, "Setting up parameters...");
        parameters.add ("dt", dt);
        parameters.add ("end_time", endTime);
        parameters.add ("dt_name", dtName);
        parameters.add ("previous_solution_name", previousSolutionName);
        parameters.add ("store_interval", storeInterval);
        parameters.add ("clone_method", "shallow_clone");
        
        dolfin::Parameters dtCoefficientTypesParameter ("dt_coefficient_types");
        for (auto& i : dtCoefficientTypes)
        {
            dtCoefficientTypesParameter.add<std::string> (i);
        }
        parameters.add (dtCoefficientTypesParameter);
        
        dolfin::Parameters previousSolutionCoefficientTypesParameter ("previous_solution_coefficient_types");
        for (auto& i : previousSolutionCoefficientTypes)
        {
            previousSolutionCoefficientTypesParameter.add<std::string> (i);
        }
        parameters.add (previousSolutionCoefficientTypesParameter);
        
        dolfin::end ();
        
        dolfin::log (dolfin::DBG, "TimeDependentDifferentialProblem object created");
    }



    TimeDependentDifferentialProblem::TimeDependentDifferentialProblem 
            (const dolfin::Mesh& mesh, 
             const dolfin::FunctionSpace& functionSpace,
             const double& dt,
             const double& endTime,
             const std::vector<std::string>& dtCoefficientTypes,
             const std::vector<std::string>& previousSolutionCoefficientTypes,
             const int& storeInterval,
             const std::shared_ptr<dcp::AbstractDifferentialProblem> timeSteppingProblem,
             const std::string& dtName,
             const std::string& previousSolutionName) :
        AbstractDifferentialProblem (mesh, functionSpace),
        timeSteppingProblem_ (timeSteppingProblem),
        solutions_ ()
    { 
        dolfin::begin (dolfin::DBG, "Building TimeDependentDifferentialProblem...");
        
        dolfin::log (dolfin::DBG, "Setting up parameters...");
        parameters.add ("dt", dt);
        parameters.add ("end_time", endTime);
        parameters.add ("dt_name", dtName);
        parameters.add ("previous_solution_name", previousSolutionName);
        parameters.add ("store_interval", storeInterval);
        parameters.add ("clone_method", "shallow_clone");
        
        dolfin::Parameters dtCoefficientTypesParameter ("dt_coefficient_types");
        for (auto& i : dtCoefficientTypes)
        {
            dtCoefficientTypesParameter.add<std::string> (i);
        }
        parameters.add (dtCoefficientTypesParameter);
        
        dolfin::Parameters previousSolutionCoefficientTypesParameter ("previous_solution_coefficient_types");
        for (auto& i : previousSolutionCoefficientTypes)
        {
            previousSolutionCoefficientTypesParameter.add<std::string> (i);
        }
        parameters.add (previousSolutionCoefficientTypesParameter);
        
        dolfin::end ();
        
        dolfin::log (dolfin::DBG, "TimeDependentDifferentialProblem object created");
    }



    TimeDependentDifferentialProblem::TimeDependentDifferentialProblem 
            (dolfin::Mesh&& mesh, 
             dolfin::FunctionSpace&& functionSpace,
             const double& dt,
             const double& endTime,
             const std::vector<std::string>& dtCoefficientTypes,
             const std::vector<std::string>& previousSolutionCoefficientTypes,
             const int& storeInterval,
             const std::shared_ptr<dcp::AbstractDifferentialProblem> timeSteppingProblem,
             const std::string& dtName,
             const std::string& previousSolutionName) :
        AbstractDifferentialProblem (mesh, functionSpace),
        timeSteppingProblem_ (timeSteppingProblem),
        solutions_ ()
    { 
        dolfin::begin (dolfin::DBG, "Building TimeDependentDifferentialProblem...");
        
        dolfin::log (dolfin::DBG, "Setting up parameters...");
        parameters.add ("dt", dt);
        parameters.add ("end_time", endTime);
        parameters.add ("dt_name", dtName);
        parameters.add ("previous_solution_name", previousSolutionName);
        parameters.add ("store_interval", storeInterval);
        parameters.add ("clone_method", "shallow_clone");
        
        dolfin::Parameters dtCoefficientTypesParameter ("dt_coefficient_types");
        for (auto& i : dtCoefficientTypes)
        {
            dtCoefficientTypesParameter.add<std::string> (i);
        }
        parameters.add (dtCoefficientTypesParameter);
        
        dolfin::Parameters previousSolutionCoefficientTypesParameter ("previous_solution_coefficient_types");
        for (auto& i : previousSolutionCoefficientTypes)
        {
            previousSolutionCoefficientTypesParameter.add<std::string> (i);
        }
        parameters.add (previousSolutionCoefficientTypesParameter);
        
        dolfin::end ();
        
        dolfin::log (dolfin::DBG, "TimeDependentDifferentialProblem object created");
    }



    /******************* GETTERS *******************/

    dcp::AbstractDifferentialProblem& TimeDependentDifferentialProblem::timeSteppingProblem ()
    {
        return *timeSteppingProblem_;
    }



    /******************* SETTERS *******************/
    void TimeDependentDifferentialProblem::
    setTimeSteppingProblem (const std::shared_ptr<dcp::AbstractDifferentialProblem> timeSteppingProblem)
    {
        timeSteppingProblem_ = timeSteppingProblem;
    }

    

    void TimeDependentDifferentialProblem::
    setInitialSolution (const dolfin::Function& initialSolution)
    {
        solution_ = initialSolution;
    }
   


    void TimeDependentDifferentialProblem::
    setCoefficient (const std::string& coefficientType, 
                    const std::shared_ptr<const dolfin::GenericFunction> coefficientValue,
                    const std::string& coefficientName)
    {
        timeSteppingProblem_->setCoefficient (coefficientType, coefficientValue, coefficientName);
    }



    void TimeDependentDifferentialProblem::
    setCoefficient (const std::string& coefficientType,
                    const std::shared_ptr<const dolfin::GenericFunction> coefficientValue,
                    const std::size_t& coefficientNumber)
    {
        timeSteppingProblem_->setCoefficient (coefficientType, coefficientValue, coefficientNumber);
    }



    void TimeDependentDifferentialProblem::
    setIntegrationSubdomains (const std::string& formType,
                              std::shared_ptr<const dolfin::MeshFunction<std::size_t>> meshFunction,
                              const dcp::SubdomainType& subdomainType)
    {
        timeSteppingProblem_->setIntegrationSubdomains (formType, meshFunction, subdomainType);
    }



    bool TimeDependentDifferentialProblem::
    addDirichletBC (const dolfin::DirichletBC& dirichletCondition, std::string bcName)
    {
        return timeSteppingProblem_->addDirichletBC (dirichletCondition, bcName);
    }



    bool TimeDependentDifferentialProblem::
    addDirichletBC (dolfin::DirichletBC&& dirichletCondition, std::string bcName)
    {
        return timeSteppingProblem_->addDirichletBC (dirichletCondition, bcName);
    }



    bool TimeDependentDifferentialProblem::
    removeDirichletBC (const std::string& bcName)
    {
        return timeSteppingProblem_->removeDirichletBC (bcName);
    }



    void TimeDependentDifferentialProblem::
    update ()
    {
        timeSteppingProblem_->update ();
    }



    /******************* METHODS *******************/
    
    void TimeDependentDifferentialProblem::
    solve () 
    {
        // ---- Problem settings ---- //
        dolfin::begin (dolfin::DBG, "Setting up time dependent problem...");
        
        dolfin::Constant dt (parameters["dt"]);
        double endTime = parameters ["end_time"];
        std::string dtName = parameters ["dt_name"];
        std::string previousSolutionName = parameters ["previous_solution_name"];
        std::vector<std::string> dtCoefficientTypes;
        parameters ("dt_coefficient_types").get_parameter_keys (dtCoefficientTypes);
        std::vector<std::string> previousSolutionCoefficientTypes;
        parameters ("previous_solution_coefficient_types").get_parameter_keys (previousSolutionCoefficientTypes);
        int storeInterval = parameters ["store_interval"];
        
        for (auto& i : dtCoefficientTypes)
        {
            timeSteppingProblem_->setCoefficient (i, dolfin::reference_to_no_delete_pointer (dt), dtName);
        } 
        
        solutions_.clear ();
        
        dolfin::end ();
        
        
        // ---- Problem solution ---- //
        dolfin::begin (dolfin::DBG, "Solving time dependent problem...");
        
        double t = 0;
        int timeStep = 0;
        
        while (t < endTime + DOLFIN_EPS)
        {
            timeStep++;
            
            dolfin::log (dolfin::INFO, "==========================");
            dolfin::log (dolfin::INFO, "Time = %f s, timestep %d", t, timeStep);
            dolfin::begin ("==========================");
            
            for (auto& i : previousSolutionCoefficientTypes)
            {
                timeSteppingProblem_->setCoefficient (i, 
                                                      dolfin::reference_to_no_delete_pointer (solution_), 
                                                      previousSolutionName);
            } 
            
            dolfin::begin ("Solving time stepping problem...");
            timeSteppingProblem_->solve ();
            dolfin::end ();
            
            solution_ = timeSteppingProblem_->solution ();
            
            if (timeStep % storeInterval == 0 )
            {
                solutions_.push_back (solution_);
            }
            
            t += dt;
            
            dolfin::end ();
        }
        
        dolfin::end ();
    }



    dcp::TimeDependentDifferentialProblem* TimeDependentDifferentialProblem::
    clone () const
    {
        dolfin::begin (dolfin::DBG, "Cloning object...");
        
        std::string cloneMethod = parameters ["clone_method"];
        
        dolfin::log (dolfin::DBG, "Clone method: %s", cloneMethod.c_str ());
        dolfin::log (dolfin::DBG, "Creating new object of type TimeDependentDifferentialProblem...");
        
        // create new object
        dcp::TimeDependentDifferentialProblem* clonedProblem = nullptr;
        if (cloneMethod == "shallow_clone")
        {
            std::vector<std::string> dtCoefficientTypes;
            parameters ("dt_coefficient_types").get_parameter_keys (dtCoefficientTypes);
            std::vector<std::string> previousSolutionCoefficientTypes;
            parameters ("previous_solution_coefficient_types").get_parameter_keys (previousSolutionCoefficientTypes);
            clonedProblem = 
                new dcp::TimeDependentDifferentialProblem (this->mesh_,
                                                           this->functionSpace_,
                                                           this->parameters ["dt"],
                                                           this->parameters ["end_time"],
                                                           dtCoefficientTypes,
                                                           previousSolutionCoefficientTypes,
                                                           this->parameters ["store_interval"],
                                                           this->timeSteppingProblem_,
                                                           this->parameters ["dt_name"],
                                                           this->parameters ["previous_solution_name"]);
        }
        else if (cloneMethod == "deep_clone")
        {
            std::vector<std::string> dtCoefficientTypes;
            parameters ("dt_coefficient_types").get_parameter_keys (dtCoefficientTypes);
            std::vector<std::string> previousSolutionCoefficientTypes;
            parameters ("previous_solution_coefficient_types").get_parameter_keys (previousSolutionCoefficientTypes);
            clonedProblem = 
                new dcp::TimeDependentDifferentialProblem (*(this->mesh_),
                                                           *(this->functionSpace_),
                                                           this->parameters ["dt"],
                                                           this->parameters ["end_time"],
                                                           dtCoefficientTypes,
                                                           previousSolutionCoefficientTypes,
                                                           this->parameters ["store_interval"],
                                                           this->timeSteppingProblem_,
                                                           this->parameters ["dt_name"],
                                                           this->parameters ["previous_solution_name"]);
        }
        else
        {
            dolfin::error ("Cannot clone linear differential problem. Unknown clone method: \"%s\"",
                           cloneMethod.c_str ());
        }
        
        //copy dirichlet boundary conditions
        dolfin::log (dolfin::DBG, "Copying Dirichlet boundary conditions...");
        for (auto &i : this->dirichletBCs_)
        {
            clonedProblem->addDirichletBC (i.second, i.first);
        }
        
        // clear parameters set of newly created object so that it can be populated by the parameters of the object
        // being created. 
        dolfin::log (dolfin::DBG, "Copying parameters to new object...");
        clonedProblem->parameters.clear ();
        clonedProblem->parameters = this->parameters;
        
        // copy solution
        dolfin::log (dolfin::DBG, "Copying solution...");
        clonedProblem->solution_ = this->solution_;
        clonedProblem->solutions_ = this->solutions_;
        
        dolfin::end ();
        
        return clonedProblem;
    }
}
