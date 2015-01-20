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

#include "IvanMovingTimeDependentProblem.h"
#include "geometry.h"
#include <dolfin/log/dolfin_log.h>
#include <time.h>

namespace Ivan
{
    /******************* CONSTRUCTORS *******************/
    MovingTimeDependentProblem::MovingTimeDependentProblem 
            (const std::shared_ptr<dolfin::Mesh> mesh, 
             const std::shared_ptr<dolfin::FunctionSpace> functionSpace,
             const double& startTime,
             const double& dt,
             const double& endTime,
             const std::vector<std::string>& dtCoefficientTypes,
             const std::vector<std::string>& previousSolutionCoefficientTypes,
             const std::shared_ptr<dcp::AbstractProblem> timeSteppingProblem,
             const int& storeInterval,
             const int& plotInterval,
             const std::string& dtName,
             const std::string& previousSolutionName) : 
        dcp::AbstractProblem (mesh, functionSpace),
        timeSteppingProblem_ (timeSteppingProblem),
        solutions_ (),
        meshManager_ (std::shared_ptr<dolfin::ALE>(new dolfin::ALE()),functionSpace),
//TODO :        displacement_ (std::shared_ptr<dolfin::GenericFunction> (new geometry::MapTgamma()))
        displacement_ (std::shared_ptr<geometry::MapTgamma> (new geometry::MapTgamma()))
    { 
        dolfin::begin (dolfin::DBG, "Building MovingTimeDependentProblem...");
        
        dolfin::log (dolfin::DBG, "Setting up parameters...");
        parameters.add ("start_time", startTime);
        parameters.add ("dt", dt);
        parameters.add ("end_time", endTime);
        parameters.add ("dt_name", dtName);
        parameters.add ("previous_solution_name", previousSolutionName);
        parameters.add ("store_interval", storeInterval);
        parameters.add ("plot_interval", plotInterval);
        parameters.add ("pause", false);
        parameters.add ("time_stepping_solution_component", -1);
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
        
        dolfin::log (dolfin::DBG, "MovingTimeDependentProblem object created");
    }



    MovingTimeDependentProblem::MovingTimeDependentProblem 
            (const dolfin::Mesh& mesh, 
             const dolfin::FunctionSpace& functionSpace,
             const double& startTime,
             const double& dt,
             const double& endTime,
             const std::vector<std::string>& dtCoefficientTypes,
             const std::vector<std::string>& previousSolutionCoefficientTypes,
             const std::shared_ptr<dcp::AbstractProblem> timeSteppingProblem,
             const int& storeInterval,
             const int& plotInterval,
             const std::string& dtName,
             const std::string& previousSolutionName) : 
        dcp::AbstractProblem (mesh, functionSpace),
        timeSteppingProblem_ (timeSteppingProblem),
        solutions_ (),
        meshManager_ (std::shared_ptr<dolfin::ALE>(new dolfin::ALE()),std::shared_ptr<dolfin::FunctionSpace>(new dolfin::FunctionSpace(functionSpace))),
//        displacement_ (std::shared_ptr<dolfin::GenericFunction> (new geometry::MapTgamma()))
        displacement_ (std::shared_ptr<geometry::MapTgamma> (new geometry::MapTgamma()))
    { 
        dolfin::begin (dolfin::DBG, "Building MovingTimeDependentProblem...");
        
        dolfin::log (dolfin::DBG, "Setting up parameters...");
        parameters.add ("start_time", startTime);
        parameters.add ("dt", dt);
        parameters.add ("end_time", endTime);
        parameters.add ("dt_name", dtName);
        parameters.add ("previous_solution_name", previousSolutionName);
        parameters.add ("store_interval", storeInterval);
        parameters.add ("plot_interval", plotInterval);
        parameters.add ("pause", false);
        parameters.add ("time_stepping_solution_component", -1);
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
        
        dolfin::log (dolfin::DBG, "MovingTimeDependentProblem object created");
    }



    MovingTimeDependentProblem::MovingTimeDependentProblem 
            (dolfin::Mesh&& mesh, 
             dolfin::FunctionSpace&& functionSpace,
             const double& startTime,
             const double& dt,
             const double& endTime,
             const std::vector<std::string>& dtCoefficientTypes,
             const std::vector<std::string>& previousSolutionCoefficientTypes,
             const std::shared_ptr<dcp::AbstractProblem> timeSteppingProblem,
             const int& storeInterval,
             const int& plotInterval,
             const std::string& dtName,
             const std::string& previousSolutionName) : 
        dcp::AbstractProblem (mesh, functionSpace),
        timeSteppingProblem_ (timeSteppingProblem),
        solutions_ (),
        meshManager_ (std::shared_ptr<dolfin::ALE>(new dolfin::ALE()),std::shared_ptr<dolfin::FunctionSpace>(&functionSpace)),
//        displacement_ (std::shared_ptr<dolfin::GenericFunction> (new geometry::MapTgamma()))
        displacement_ (std::shared_ptr<geometry::MapTgamma> (new geometry::MapTgamma()))
    { 
        dolfin::begin (dolfin::DBG, "Building MovingTimeDependentProblem...");
        
        dolfin::log (dolfin::DBG, "Setting up parameters...");
        parameters.add ("start_time", startTime);
        parameters.add ("dt", dt);
        parameters.add ("end_time", endTime);
        parameters.add ("dt_name", dtName);
        parameters.add ("previous_solution_name", previousSolutionName);
        parameters.add ("store_interval", storeInterval);
        parameters.add ("plot_interval", plotInterval);
        parameters.add ("pause", false);
        parameters.add ("time_stepping_solution_component", -1);
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
        
        dolfin::log (dolfin::DBG, "MovingTimeDependentProblem object created");
    }



    /******************* GETTERS *******************/
    dcp::AbstractProblem& MovingTimeDependentProblem::timeSteppingProblem ()
    {
        return *timeSteppingProblem_;
    }

    

    const dolfin::Function& MovingTimeDependentProblem::solution () const
    {
        return solution_;
    }

    

    const std::vector<dolfin::Function>& MovingTimeDependentProblem::solutionsVector () const
    {
        return solutions_;
    }



    /******************* SETTERS *******************/
    void MovingTimeDependentProblem::
    setTimeSteppingProblem (const std::shared_ptr<dcp::AbstractProblem> timeSteppingProblem)
    {
        timeSteppingProblem_ = timeSteppingProblem;
    }

    

    void MovingTimeDependentProblem::
    setInitialSolution (const dolfin::Function& initialSolution)
    {
        solution_ = initialSolution;
    }

    

    void MovingTimeDependentProblem::
    setInitialSolution (const dolfin::Expression& initialSolution)
    {
        solution_ = initialSolution;
    }
   


    void MovingTimeDependentProblem::
    setCoefficient (const std::string& coefficientType, 
                    const std::shared_ptr<const dolfin::GenericFunction> coefficientValue,
                    const std::string& coefficientName)
    {
        timeSteppingProblem_->setCoefficient (coefficientType, coefficientValue, coefficientName);
    }



    void MovingTimeDependentProblem::
    setCoefficient (const std::string& coefficientType,
                    const std::shared_ptr<const dolfin::GenericFunction> coefficientValue,
                    const std::size_t& coefficientNumber)
    {
        timeSteppingProblem_->setCoefficient (coefficientType, coefficientValue, coefficientNumber);
    }



    void MovingTimeDependentProblem::
    setIntegrationSubdomains (const std::string& formType,
                              std::shared_ptr<const dolfin::MeshFunction<std::size_t>> meshFunction,
                              const dcp::SubdomainType& subdomainType)
    {
        timeSteppingProblem_->setIntegrationSubdomains (formType, meshFunction, subdomainType);
    }



    bool MovingTimeDependentProblem::
    addDirichletBC (const dolfin::DirichletBC& dirichletCondition, std::string bcName)
    {
        return timeSteppingProblem_->addDirichletBC (dirichletCondition, bcName);
    }



    bool MovingTimeDependentProblem::
    addDirichletBC (dolfin::DirichletBC&& dirichletCondition, std::string bcName)
    {
        return timeSteppingProblem_->addDirichletBC (dirichletCondition, bcName);
    }



    bool MovingTimeDependentProblem::
    removeDirichletBC (const std::string& bcName)
    {
        return timeSteppingProblem_->removeDirichletBC (bcName);
    }



    void MovingTimeDependentProblem::
    update ()
    {
        timeSteppingProblem_->update ();
    }



    /******************* METHODS *******************/
    void Ivan::MovingTimeDependentProblem::
    solve () 
    {
        // ---- Problem settings ---- //
        dolfin::begin (dolfin::DBG, "Setting up time dependent problem...");
        
        // get parameters' values
        dolfin::Constant dt (parameters["dt"]);
        double endTime = parameters ["end_time"];
        std::string dtName = parameters ["dt_name"];
        std::string previousSolutionName = parameters ["previous_solution_name"];
        std::vector<std::string> dtCoefficientTypes;
        parameters ("dt_coefficient_types").get_parameter_keys (dtCoefficientTypes);
        std::vector<std::string> previousSolutionCoefficientTypes;
        parameters ("previous_solution_coefficient_types").get_parameter_keys (previousSolutionCoefficientTypes);
        int storeInterval = parameters ["store_interval"];
        int plotInterval = parameters ["plot_interval"];
        bool pause = parameters ["pause"];
        int timeSteppingSolutionComponent = parameters ["time_stepping_solution_component"];
        
        for (auto& i : dtCoefficientTypes)
        {
            timeSteppingProblem_->setCoefficient (i, dolfin::reference_to_no_delete_pointer (dt), dtName);
        } 
        
        solutions_.clear ();
        
        dolfin::end ();
        
        
        // ---- Problem solution ---- //
        dolfin::begin (dolfin::DBG, "Solving time dependent problem...");
        
        double t = parameters ["start_time"];
        int timeStep = 0;
        displacement_->initTime(t);
        
        // save initial solution in member vector
        solutions_.push_back (solution_);
        
        while (t < endTime + DOLFIN_EPS)
        {
            timeStep++;
            t += dt;
            
            dolfin::begin (dolfin::INFO, "===== Time = %f s, timestep %d =====", t, timeStep);
            
            for (auto& i : previousSolutionCoefficientTypes)
            {
                if (timeSteppingSolutionComponent >= 0)
                {
                    timeSteppingProblem_->setCoefficient 
                        (i, 
                         dolfin::reference_to_no_delete_pointer (solution_ [timeSteppingSolutionComponent]), 
                         previousSolutionName);
                }
                else
                {
                    timeSteppingProblem_->setCoefficient 
                        (i, 
                         dolfin::reference_to_no_delete_pointer (solution_), 
                         previousSolutionName);
                }
            } 
            
            dolfin::begin (dolfin::INFO, "Solving time stepping problem...");
            timeSteppingProblem_->solve ();
            dolfin::end ();
            
            solution_ = timeSteppingProblem_->solution ();
            
            if (storeInterval > 0 && timeStep % storeInterval == 0)
            {
                dolfin::log (dolfin::DBG, "Saving time stepping problem solution in solutions vector...");
                solutions_.push_back (solution_);
            }
            
            if (plotInterval > 0 && timeStep % plotInterval == 0)
            {
                dolfin::log (dolfin::DBG, "Plotting time stepping problem solution...");
                dolfin::plot (solution_);
                
                if (pause)
                {
                    dolfin::interactive ();
                }
            }

            displacement_->setTime(t);
            meshManager_.moveMesh(*displacement_);
            
            dolfin::end ();

            sleep(1);
        }
        
        // save final solution in member vector if it has not been saved already
        if (! (storeInterval > 0 && timeStep % storeInterval == 0) )
        {
            dolfin::log (dolfin::DBG, "Saving time stepping problem solution in solutions vector...");
            solutions_.push_back (solution_);
        }
        
        dolfin::end ();
    }



    MovingTimeDependentProblem* MovingTimeDependentProblem::
    clone () const
    {
        dolfin::begin (dolfin::DBG, "Cloning object...");
        
        std::string cloneMethod = parameters ["clone_method"];
        
        dolfin::log (dolfin::DBG, "Clone method: %s", cloneMethod.c_str ());
        dolfin::log (dolfin::DBG, "Creating new object of type MovingTimeDependentProblem...");
        
        // create new object
        MovingTimeDependentProblem* clonedProblem = nullptr;
        if (cloneMethod == "shallow_clone")
        {
            std::vector<std::string> dtCoefficientTypes;
            parameters ("dt_coefficient_types").get_parameter_keys (dtCoefficientTypes);
            std::vector<std::string> previousSolutionCoefficientTypes;
            parameters ("previous_solution_coefficient_types").get_parameter_keys (previousSolutionCoefficientTypes);
            clonedProblem = 
                new MovingTimeDependentProblem (this->mesh_,
                                                           this->functionSpace_,
                                                           this->parameters ["start_time"],
                                                           this->parameters ["dt"],
                                                           this->parameters ["end_time"],
                                                           dtCoefficientTypes,
                                                           previousSolutionCoefficientTypes,
                                                           this->timeSteppingProblem_,
                                                           this->parameters ["store_interval"],
                                                           this->parameters ["plot_interval"],
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
                new MovingTimeDependentProblem (*(this->mesh_),
                                                           *(this->functionSpace_),
                                                           this->parameters ["start_time"],
                                                           this->parameters ["dt"],
                                                           this->parameters ["end_time"],
                                                           dtCoefficientTypes,
                                                           previousSolutionCoefficientTypes,
                                                           this->timeSteppingProblem_,
                                                           this->parameters ["store_interval"],
                                                           this->parameters ["plot_interval"],
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
