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

#include <differential_problems/TimeDependentProblem.h>
#include <dolfin/log/dolfin_log.h>
#include <regex>

namespace dcp
{
    /******************* CONSTRUCTORS *******************/
    TimeDependentProblem::TimeDependentProblem 
        (const std::shared_ptr<dcp::AbstractProblem> timeSteppingProblem,
         const double& startTime,
         const double& dt,
         const double& endTime,
         std::initializer_list<std::string> dtCoefficientTypes,
         std::initializer_list<std::string> previousSolutionCoefficientTypes,
         const int& storeInterval,
         const int& plotInterval,
         const std::string& dtName,
         const std::string& previousSolutionName) 
        : 
            AbstractProblem (timeSteppingProblem->functionSpace ()),
            timeSteppingProblem_ (timeSteppingProblem),
            t_ (startTime),
            startTime_ (startTime),
            dt_ (dt),
            endTime_ (endTime)
    { 
        dolfin::begin (dolfin::DBG, "Building TimeDependentProblem...");
        
        solution_.emplace_back (dolfin::Function (timeSteppingProblem->functionSpace ()));
            
        dolfin::log (dolfin::DBG, "Setting up parameters...");
        parameters.add ("problem_type", "time_dependent");
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
        
        dolfin::log (dolfin::DBG, "TimeDependentProblem object created");
    }



    /******************* GETTERS *******************/
    std::shared_ptr<const dolfin::Mesh> TimeDependentProblem::mesh () const
    {
        return timeSteppingProblem_ -> mesh ();
    }



    std::shared_ptr<dolfin::FunctionSpace> TimeDependentProblem::functionSpace () const
    {
        return timeSteppingProblem_ -> functionSpace ();
    }



    const dolfin::DirichletBC& TimeDependentProblem::dirichletBC (const std::string& bcName) const
    {
        return timeSteppingProblem_ -> dirichletBC (bcName);
    }

    

    const std::map<std::string, dolfin::DirichletBC>& TimeDependentProblem::dirichletBCs () const
    {
        return timeSteppingProblem_ -> dirichletBCs ();
    }
    


    const dolfin::Function& TimeDependentProblem::solution () const
    {
        return solution_.back ();
    }



    const std::vector<dolfin::Function>& TimeDependentProblem::solutionsVector () const
    {
        return solution_;
    }
    


    const double& TimeDependentProblem::time () const
    {
        return t_;
    }
    


    double& TimeDependentProblem::startTime ()
    {
        return startTime_;
    }
    


    double& TimeDependentProblem::dt ()
    {
        return dt_;
    }



    double& TimeDependentProblem::endTime ()
    {
        return endTime_;
    }
    


    dcp::AbstractProblem& TimeDependentProblem::timeSteppingProblem ()
    {
        return *timeSteppingProblem_;
    }



    /******************* SETTERS *******************/
    void TimeDependentProblem::setInitialSolution (const dolfin::Function& initialSolution)
    {
        solution_.back () = initialSolution;
    }

    

    void TimeDependentProblem::setInitialSolution (const dolfin::Expression& initialSolution)
    {
        solution_.back () = initialSolution;
    }
   


    void TimeDependentProblem::setCoefficient (const std::string& coefficientType, 
                                               const std::shared_ptr<const dolfin::GenericFunction> coefficientValue,
                                               const std::string& coefficientName)
    {
        timeSteppingProblem_->setCoefficient (coefficientType, coefficientValue, coefficientName);
    }



    void TimeDependentProblem::setCoefficient (const std::string& coefficientType,
                                               const std::shared_ptr<const dolfin::GenericFunction> coefficientValue,
                                               const std::size_t& coefficientNumber)
    {
        timeSteppingProblem_->setCoefficient (coefficientType, coefficientValue, coefficientNumber);
    }



    void TimeDependentProblem::
    setIntegrationSubdomains (const std::string& formType,
                              std::shared_ptr<const dolfin::MeshFunction<std::size_t>> meshFunction,
                              const dcp::SubdomainType& subdomainType)
    {
        timeSteppingProblem_->setIntegrationSubdomains (formType, meshFunction, subdomainType);
    }



    bool TimeDependentProblem::addDirichletBC (const dolfin::GenericFunction& condition, 
                                               const dolfin::SubDomain& boundary,
                                               std::string bcName)
    {
        return timeSteppingProblem_->addDirichletBC (condition, boundary, bcName);
    }
    


    bool TimeDependentProblem::addDirichletBC (const dolfin::GenericFunction& condition, 
                                               const dolfin::SubDomain& boundary, 
                                               const std::size_t& component,
                                               std::string bcName)
    {
        return timeSteppingProblem_->addDirichletBC (condition, boundary, component, bcName);
    }
    


    bool TimeDependentProblem::addDirichletBC (std::shared_ptr<const dolfin::GenericFunction> condition, 
                                               std::shared_ptr<const dolfin::SubDomain> boundary,
                                               std::string bcName)
    {
        return timeSteppingProblem_->addDirichletBC (condition, boundary, bcName);
    }
    


    bool TimeDependentProblem::addDirichletBC (std::shared_ptr<const dolfin::GenericFunction> condition, 
                                               std::shared_ptr<const dolfin::SubDomain> boundary, 
                                               const std::size_t& component,
                                               std::string bcName)
    {
        return timeSteppingProblem_->addDirichletBC (condition, boundary, component, bcName);
    }
    


    bool TimeDependentProblem::addDirichletBC (const dolfin::DirichletBC& dirichletCondition, 
                                               std::string bcName)
    {
        return timeSteppingProblem_->addDirichletBC (dirichletCondition, bcName);
    }



    bool TimeDependentProblem::addDirichletBC (dolfin::DirichletBC&& dirichletCondition, 
                                               std::string bcName)
    {
        return timeSteppingProblem_->addDirichletBC (dirichletCondition, bcName);
    }



    bool TimeDependentProblem::removeDirichletBC (const std::string& bcName)
    {
        return timeSteppingProblem_->removeDirichletBC (bcName);
    }



    void TimeDependentProblem::update ()
    {
        timeSteppingProblem_->update ();
    }



    /******************* METHODS *******************/
    bool TimeDependentProblem::isFinished ()
    {
        return (dt_ > 0) ? (t_ >= endTime_) : (t_ <= endTime_);
    }
    


    void TimeDependentProblem::clear ()
    {
        solution_.clear ();
        solution_.emplace_back (dolfin::Function (timeSteppingProblem_->functionSpace ()));
        
        t_ = startTime_;
    }
    


    void TimeDependentProblem::step ()
    {
        solve ("step");
    }
    

    
    void TimeDependentProblem::solve (const std::string& type) 
    {
        // parse solve type
        if (type != "default" && type != "step" && type != "clear_default" && type != "clear_step")
        {
            dolfin::dolfin_error ("dcp: TimeDependentProblem.h", 
                                  "solve",
                                  "Unknown solve type \"%s\" requested",
                                  type.c_str ());
        }
        
        // check if we are already at the end of time loop
        if (isFinished ())
        {
            printFinishedWarning ();
            return;
        }

        dolfin::log (dolfin::DBG, "Solve type: %s", type.c_str ());
            
        if (std::regex_match (type, std::regex (".*clear.*")))
        {
            clear ();
        }
        
        // ---- Problem settings ---- //
        dolfin::begin (dolfin::DBG, "Setting up time dependent problem...");
        
        // flag to check if the solver must perform just one time step
        bool oneStepRequested = std::regex_match (type, std::regex (".*step.*"));
        
        // get parameters' values
        dolfin::Constant dt (dt_);
        std::string dtName = parameters ["dt_name"];
        std::vector<std::string> dtCoefficientTypes;
        parameters ("dt_coefficient_types").get_parameter_keys (dtCoefficientTypes);
        int storeInterval = parameters ["store_interval"];
        int plotInterval = parameters ["plot_interval"];
        bool pause = parameters ["pause"];
        
        for (auto& i : dtCoefficientTypes)
        {
            timeSteppingProblem_->setCoefficient (i, dolfin::reference_to_no_delete_pointer (dt), dtName);
        } 
        
        dolfin::end ();
        
        
        // ---- Problem solution ---- //
        dolfin::begin (dolfin::DBG, "Solving time dependent problem...");
        
        
        // function used to step through the time loop
        dolfin::Function tmpSolution = solution_.back ();
        
        // plot initial solution
        dolfin::plot (tmpSolution, "Initial solution");
        if (pause)
        {
            dolfin::interactive ();
        }
        
        // start time loop. The loop variable oneTimeStep is true only if the solve type requested is "step" and
        // the first step has already been peformed
        int timeStep = 0;
        bool timeStepFlag = false;
        while (!isFinished () && timeStepFlag == false)
        {
            timeStep++;
            
            dolfin::begin (dolfin::INFO, "===== Timestep %d =====", timeStep);
            
            advanceTime (tmpSolution);
            
            dolfin::begin (dolfin::INFO, "Solving time stepping problem...");
            timeSteppingProblem_->solve ();
            dolfin::end ();
            
            tmpSolution = timeSteppingProblem_->solution ();
            
            // save solution in solution_ according to time step and store interval.
            if (storeInterval > 0 && timeStep % storeInterval == 0)
            {
                dolfin::log (dolfin::DBG, "Saving time stepping problem solution in solutions vector...");
                solution_.push_back (tmpSolution);
            }
            
            if (plotInterval > 0 && timeStep % plotInterval == 0)
            {
                dolfin::log (dolfin::DBG, "Plotting time stepping problem solution...");
                dolfin::plot (tmpSolution, "Time = " + std::to_string (t_));
                
                if (pause)
                {
                    dolfin::interactive ();
                }
            }
            
            if (oneStepRequested == true)
            {
                timeStepFlag = true;
            }
            
            dolfin::end ();
        }
        
        // At this point, we just need to make sure that the solution on the last iteration was saved even though
        // timeStep % storeInterval != 0 (but it must not be saved twice!)
        if (!(storeInterval > 0 && timeStep % storeInterval == 0))
        {
            dolfin::log (dolfin::DBG, "Saving last time step solution in solutions vector...");
            solution_.push_back (tmpSolution);
        }
        
        dolfin::end ();
    }



    dcp::TimeDependentProblem* TimeDependentProblem::
    clone () const
    {
        dolfin::begin (dolfin::DBG, "Cloning object...");
        
        std::string cloneMethod = parameters ["clone_method"];
        
        dolfin::log (dolfin::DBG, "Clone method: %s", cloneMethod.c_str ());
        dolfin::log (dolfin::DBG, "Creating new object of type TimeDependentProblem...");
        
        // create new object
        dcp::TimeDependentProblem* clonedProblem = nullptr;
        if (cloneMethod == "shallow_clone")
        {
            // note that we pass an empty initializer_list to the constructor as dtCoefficientTypes and 
            // previousSolutionCoefficientTypes, because they will be copied when the parameters are copied anyway
            clonedProblem = 
                new dcp::TimeDependentProblem (this->timeSteppingProblem_,
                                               startTime_,
                                               dt_,
                                               endTime_,
                                               {},
                                               {},
                                               this->parameters ["store_interval"],
                                               this->parameters ["plot_interval"],
                                               this->parameters ["dt_name"],
                                               this->parameters ["previous_solution_name"]);
        }
        else if (cloneMethod == "deep_clone")
        {
            // note that we pass an empty initializer_list to the constructor as dtCoefficientTypes and 
            // previousSolutionCoefficientTypes, because they will be copied when the parameters are copied anyway
            clonedProblem = 
                new dcp::TimeDependentProblem (this->timeSteppingProblem_,
                                               startTime_,
                                               dt_,
                                               endTime_,
                                               {},
                                               {},
                                               this->parameters ["store_interval"],
                                               this->parameters ["plot_interval"],
                                               this->parameters ["dt_name"],
                                               this->parameters ["previous_solution_name"]);
        }
        else
        {
            dolfin::dolfin_error ("dcp: TimeDependentProblem.cpp",
                                  "clone",
                                  "Cannot clone time dependent differential problem. Unknown clone method: \"%s\"",
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
        
        dolfin::end ();
        
        return clonedProblem;
    }
    


    /******************* PROTECTED METHODS *******************/
    void TimeDependentProblem::advanceTime (const dolfin::Function& previousSolution)
    {
        t_ += dt_;
        
        dolfin::log (dolfin::INFO, "TIME = %f s", t_);

        std::string previousSolutionName = parameters ["previous_solution_name"];
        std::vector<std::string> previousSolutionCoefficientTypes;
        parameters ("previous_solution_coefficient_types").get_parameter_keys (previousSolutionCoefficientTypes);
        int timeSteppingSolutionComponent = parameters ["time_stepping_solution_component"];
        
        dolfin::log (dolfin::DBG, "Setting previous solution coefficients...");
        for (auto& i : previousSolutionCoefficientTypes)
        {
            if (timeSteppingSolutionComponent >= 0)
            {
                timeSteppingProblem_->setCoefficient 
                    (i, 
                     dolfin::reference_to_no_delete_pointer (previousSolution [timeSteppingSolutionComponent]), 
                     previousSolutionName);
            }
            else
            {
                timeSteppingProblem_->setCoefficient 
                    (i, 
                     dolfin::reference_to_no_delete_pointer (previousSolution), 
                     previousSolutionName);
            }
        } 
    }
    


    void TimeDependentProblem::printFinishedWarning ()
    {
        dolfin::warning ("No time iteration performed in solve() function. End time already reached.");
    }
}
