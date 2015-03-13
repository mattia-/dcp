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

#include <dcp/differential_problems/TimeDependentProblem.h>
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
         const unsigned int& nTimeSchemeSteps)
        : 
            AbstractProblem (timeSteppingProblem->functionSpace ()),
            timeSteppingProblem_ (timeSteppingProblem),
            timeDependentCoefficients_ (),
            timeDependentDirichletBCs_ (),
            t_ (startTime),
            startTime_ (startTime),
            dt_ (dt),
            endTime_ (endTime),
            nTimeSchemeSteps_ (nTimeSchemeSteps),
            timeDependentDirichletBCsCounter_ (0)
    { 
        dolfin::begin (dolfin::DBG, "Building TimeDependentProblem...");
        
        // we need a solution for every step in the time scheme! 
        // And the time should be t0, t0+dt, t0+2dt ...
        dolfin::begin (dolfin::DBG, "Creating initial solutions...");
        unsigned int stepNumber;
        for (stepNumber = 0; stepNumber < nTimeSchemeSteps_; ++stepNumber)
        {
            solution_.emplace_back (std::make_pair (t_ + stepNumber * dt, 
                                                    dolfin::Function (timeSteppingProblem_->functionSpace ())));
        }
        dolfin::log (dolfin::DBG, "Created %d initial solutions", stepNumber + 1);
        
        // set the correct value for t_, since it was not incremented during the previous loop.
        // Use stepNumber - 1 since t_ will be incremented at the *beginning* of the time loop, not at the end
        t_ += (stepNumber - 1) * dt;
            
        dolfin::end ();
        
        dolfin::log (dolfin::DBG, "Setting up parameters...");
        parameters.add ("problem_type", "time_dependent");
        parameters.add ("dt_name", "dt");
        parameters.add ("previous_solution_name", "u_old");
        parameters.add ("store_interval", 1);
        parameters.add ("plot_interval", 0);
        parameters.add ("time_stepping_solution_component", -1);
        parameters.add ("pause", false);
        parameters.add ("previous_solution_is_set_externally", false);
        
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
        
        parameters ["plot_title"] = "";

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
        return solution_.back ().second;
    }



    const std::vector <std::pair <double, dolfin::Function> >& TimeDependentProblem::solutionsVector () const
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
    void TimeDependentProblem::setInitialSolution (const dolfin::Function& initialSolution, 
                                                   const unsigned int& stepNumber)
    {
        if (stepNumber > nTimeSchemeSteps_)
        {
            dolfin::warning ("initial solution not set. Requested time step number is greater than problem's scheme's number of time steps");
            return;
        }
        
        auto solutionsIterator = solution_.rbegin ();
        
        // the solution we want to set is identified by stepNumber. If it is equal to 1, we must set the last element
        // in solution_, if it is equal to 2 the last but one element and so on. That's why we subtract 1 from 
        // stepNumber in the next instruction
        (solutionsIterator + (stepNumber - 1)) -> second = initialSolution;
    }

    

    void TimeDependentProblem::setInitialSolution (const dolfin::Expression& initialSolution,
                                                   const unsigned int& stepNumber)
    {
        if (stepNumber > nTimeSchemeSteps_)
        {
            dolfin::warning ("initial solution not set. Requested time step number is greater than problem's scheme's number of time steps");
            return;
        }
        
        auto solutionsIterator = solution_.rbegin ();
        
        // the solution we want to set is identified by stepNumber. If it is equal to 1, we must set the last element
        // in solution_, if it is equal to 2 the last but one element and so on. That's why we subtract 1 from 
        // stepNumber in the next instruction
        (solutionsIterator + (stepNumber - 1)) -> second = initialSolution;
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
    setIntegrationSubdomain (const std::string& formType,
                              std::shared_ptr<const dolfin::MeshFunction<std::size_t>> meshFunction,
                              const dcp::SubdomainType& subdomainType)
    {
        timeSteppingProblem_->setIntegrationSubdomain (formType, meshFunction, subdomainType);
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

    
        
    bool TimeDependentProblem::addTimeDependentDirichletBC (const dcp::TimeDependentExpression& condition, 
                                                            const dcp::Subdomain& boundary,
                                                            std::string bcName)
    {
        if (bcName.empty ())
        {
            bcName = "time_dependent_dirichlet_condition_" + std::to_string (timeDependentDirichletBCsCounter_);
            timeDependentDirichletBCsCounter_++;
        }

        dolfin::begin (dolfin::DBG, 
                       "Adding dirichlet boundary condition to time dependent boundary conditions map with name \"%s\"...",
                       bcName.c_str ());

        auto result = timeDependentDirichletBCs_.emplace 
            (bcName, 
             std::make_tuple (std::shared_ptr<dcp::TimeDependentExpression> (condition.clone ()), 
                              std::shared_ptr<dcp::Subdomain> (boundary.clone ()), 
                              -1)
             );

        bool bcWasAdded = result.second;

        if (bcWasAdded == false)
        {
            dolfin::warning ("DirichletBC object not inserted because key \"%s\" already in map",
                             bcName.c_str ());
        }
        else
        {
            dolfin::begin (dolfin::DBG, 
                           "Adding dirichlet boundary condition to time stepping problem boundary conditions map with name \"%s\"...",
                           bcName.c_str ());
            bcWasAdded = bcWasAdded && addDirichletBC (std::get<0> ((result.first)->second),
                                                       std::get<1> ((result.first)->second),
                                                       bcName);
            dolfin::end ();
        }
        
        dolfin::end ();

        return bcWasAdded;
    }



    bool TimeDependentProblem::addTimeDependentDirichletBC (const dcp::TimeDependentExpression& condition, 
                                                            const dcp::Subdomain& boundary,
                                                            const std::size_t& component,
                                                            std::string bcName)
    {
        if (bcName.empty ())
        {
            bcName = "time_dependent_dirichlet_condition_" + std::to_string (timeDependentDirichletBCsCounter_);
            timeDependentDirichletBCsCounter_++;
        }

        dolfin::begin (dolfin::DBG, 
                       "Adding dirichlet boundary condition to time dependent boundary conditions map with name \"%s\"...",
                       bcName.c_str ());

        auto result = timeDependentDirichletBCs_.emplace 
            (bcName, 
             std::make_tuple (std::shared_ptr<dcp::TimeDependentExpression> (condition.clone ()), 
                              std::shared_ptr<dcp::Subdomain> (boundary.clone ()), 
                              component)
             );

        bool bcWasAdded = result.second;

        if (bcWasAdded == false)
        {
            dolfin::warning ("DirichletBC object not inserted because key \"%s\" already in map",
                             bcName.c_str ());
        }
        else
        {
            dolfin::begin (dolfin::DBG, 
                           "Adding dirichlet boundary condition to time stepping problem boundary conditions map with name \"%s\"...",
                           bcName.c_str ());
            bcWasAdded = bcWasAdded && addDirichletBC (std::get<0> ((result.first)->second),
                                                       std::get<1> ((result.first)->second),
                                                       component,
                                                       bcName);
            dolfin::end ();
        }

        dolfin::end ();

        return result.second;
    }
    
    

    bool TimeDependentProblem::removeTimeDependentDirichletBC (const std::string& bcName)
    {
        dolfin::begin (dolfin::DBG, 
                     "Removing dirichlet boundary condition \"%s\" from time dependent boundary conditions map...", 
                     bcName.c_str ());
        std::size_t nErasedElements = timeDependentDirichletBCs_.erase (bcName);

        if (nErasedElements == 0)
        {
            dolfin::warning ("Dirichlet boundary condition \"%s\" not found in map", 
                             bcName.c_str ());
        }
        
        dolfin::begin (dolfin::DBG, 
                       "Removing dirichlet boundary condition \"%s\" from time stepping problem boundary conditions map...", 
                       bcName.c_str ());
        // call removeDirichletBC to remove the dirichlet bc also from timeSteppingProblem_
        // We add it to the existing value of nErasedElements so that we know wether both removal operations went fine
        nErasedElements += removeDirichletBC (bcName);
        dolfin::end ();

        dolfin::end ();
        
        // remember: both removal were ok if nErasedElements was equal to 1 both times!
        return nErasedElements == 2? true : false;
    }



    bool TimeDependentProblem::addTimeDependentCoefficient (const std::string& coefficientName, 
                                                            const std::string& coefficientType,
                                                            std::shared_ptr<dcp::TimeDependentExpression> expression)
    {
        dolfin::begin (dolfin::DBG, 
                       "Inserting time dependent coefficient in map with name \"%s\" and type \"%s\"...",
                       coefficientName.c_str (),
                       coefficientType.c_str ());
       
        auto coefficientID = std::make_pair (coefficientName, coefficientType);
        auto result = timeDependentCoefficients_.insert (std::make_pair (coefficientID, expression));
        
        if (result.second == false)
        {
            dolfin::warning ("Time dependent coefficient not inserted because key is already present in map");
        }
        
        dolfin::end ();
        
        return result.second;
    }



    bool TimeDependentProblem::removeTimeDependentCoefficient (const std::string& coefficientName,
                                                               const std::string& coefficientType)
    {
        dolfin::begin (dolfin::DBG, 
                       "Removing time dependent coefficient with name \"%s\" and type \"%s\" from map...",
                       coefficientName.c_str (),
                       coefficientType.c_str ());
       
        auto coefficientID = std::make_pair (coefficientName, coefficientType);
        
        std::size_t nErasedElements = timeDependentCoefficients_.erase (coefficientID);
        
        if (nErasedElements == 0)
        {
            dolfin::warning ("Time dependent coefficient not removed: key was not found in map");
        }
        
        dolfin::end ();
        
        return nErasedElements == 1? true : false;
    }
            


    void TimeDependentProblem::update ()
    {
        timeSteppingProblem_->update ();
    }



    /******************* METHODS *******************/
    bool TimeDependentProblem::isFinished ()
    {
        return (dt_ > 0) ? (t_ >= endTime_ - DOLFIN_EPS) : (t_ <= endTime_ + DOLFIN_EPS);
    }
    


    void TimeDependentProblem::clear ()
    {
        // reset t_
        t_ = startTime_;
        
        // clear solutions vector
        solution_.clear ();
        solution_.emplace_back (std::make_pair (t_, dolfin::Function (timeSteppingProblem_->functionSpace ())));
        
        // reset time dependent Dirichlet BCs
        for (auto bcIterator = timeDependentDirichletBCs_.begin (); 
             bcIterator != timeDependentDirichletBCs_.end (); 
             bcIterator++)
        {
            resetTimeDependentDirichletBC (bcIterator);
        }
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
        int plotComponent = parameters ["plot_component"];
        bool pause = parameters ["pause"];
        
        for (auto& i : dtCoefficientTypes)
        {
            timeSteppingProblem_->setCoefficient (i, dolfin::reference_to_no_delete_pointer (dt), dtName);
        } 
        
        dolfin::end ();
        
        
        // ---- Problem solution ---- //
        dolfin::begin (dolfin::DBG, "Solving time dependent problem...");
        
        
        // function used to step through the time loop
        dolfin::Function tmpSolution = solution_.back ().second;
        
        std::shared_ptr<dolfin::VTKPlotter> plotter;
        
        // start time loop. The loop variable timeStepFlag is true only if the solve type requested is "step" and
        // the first step has already been peformed
        int timeStep = 0;
        bool timeStepFlag = false;
        while (!isFinished () && timeStepFlag == false)
        {
            timeStep++;
            
            if (oneStepRequested == false)
            {
                dolfin::begin (dolfin::INFO, "===== Timestep %d =====", timeStep);
            }
            
            advanceTime ();
            
            setTimeDependentCoefficients ();
            
            setTimeDependentDirichletBCs ();
            
            dolfin::begin (dolfin::INFO, "Solving time stepping problem...");
            timeSteppingProblem_->solve ();
            dolfin::end ();
            
            tmpSolution = timeSteppingProblem_->solution ();
            
            // save solution in solution_ according to time step and store interval.
            storeSolution (tmpSolution, timeStep, storeInterval);
            
            // plot solution according to time step and plot interval
            plotSolution (tmpSolution, timeStep, plotInterval, plotComponent, pause);
            
            if (oneStepRequested == true)
            {
                timeStepFlag = true;
            }
             
            // dolfin::end matching dolfin::begin inside if statement
            if (oneStepRequested == false)
            {
                dolfin::end ();
            }
        }
        
        // At this point, we just need to make sure that the solution on the last iteration was saved even though
        // timeStep % storeInterval != 0 (but it must not be saved twice!)
        storeLastStepSolution (tmpSolution, timeStep, storeInterval);
        
        dolfin::end ();
    }
    
    

    void TimeDependentProblem::plotSolution ()
    {
        bool pause = parameters ["pause"];
        int plotComponent = parameters ["plot_component"];
        std::string plotTitle = parameters ["plot_title"];
        
        // if plotTitle is not empty, we need to prepend ", " so that the plot title is readable
        if (!plotTitle.empty ())
        {
            plotTitle = ", " + plotTitle;
        }
        
        // auxiliary variable, to enhance readability
        std::shared_ptr<dolfin::Function> functionToPlot;
        
        dolfin::begin (dolfin::DBG, "Plotting...");
        
        for (auto& timeSolutionPair : solution_)
        {
            double time = timeSolutionPair.first;
            // get right function to plot
            if (plotComponent == -1)
            {
                functionToPlot = dolfin::reference_to_no_delete_pointer (timeSolutionPair.second);
                dolfin::log (dolfin::DBG, "Plotting time stepping problem solution, all components...");
            }
            else
            {
                functionToPlot = dolfin::reference_to_no_delete_pointer (timeSolutionPair.second [plotComponent]);
                dolfin::log (dolfin::DBG, "Plotting time stepping problem solution, component %d...", plotComponent);
            }

            // actual plotting
            if (solutionPlotter_ == nullptr)
            {
                dolfin::log (dolfin::DBG, "Plotting in new dolfin::VTKPlotter object...");
                solutionPlotter_ = dolfin::plot (functionToPlot, "Time = " + std::to_string (time) + plotTitle);
            }
            else if (! solutionPlotter_ -> is_compatible (functionToPlot))
            {
                dolfin::log (dolfin::DBG, "Existing plotter is not compatible with object to be plotted.");
                dolfin::log (dolfin::DBG, "Creating new dolfin::VTKPlotter object...");
                solutionPlotter_ = dolfin::plot (functionToPlot, "Time = " + std::to_string (time) + plotTitle);
            }
            else 
            {
                solutionPlotter_ -> parameters ["title"] = std::string ("Time = " + std::to_string (time) + plotTitle);
                solutionPlotter_ -> plot (functionToPlot);
            }

            if (pause)
            {
                dolfin::interactive ();
            }
        }
        
        dolfin::end ();
    }



    dcp::TimeDependentProblem* TimeDependentProblem::clone () const
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
                new dcp::TimeDependentProblem (this->timeSteppingProblem_, startTime_, dt_, endTime_, {}, {});
            clonedProblem->timeDependentDirichletBCs_ = this->timeDependentDirichletBCs_;
        }
        else if (cloneMethod == "deep_clone")
        {
            // note that we pass an empty initializer_list to the constructor as dtCoefficientTypes and 
            // previousSolutionCoefficientTypes, because they will be copied when the parameters are copied anyway
            clonedProblem = 
                new dcp::TimeDependentProblem (this->timeSteppingProblem_, startTime_, dt_, endTime_, {}, {});
        }
        else
        {
            dolfin::dolfin_error ("dcp: TimeDependentProblem.cpp",
                                  "clone",
                                  "Cannot clone time dependent differential problem. Unknown clone method: \"%s\"",
                                  cloneMethod.c_str ());
            for (auto& bc : this->timeDependentDirichletBCs_)
            {
                clonedProblem->addTimeDependentDirichletBC (*(std::get<0> (bc.second)),
                                                            *(std::get<1> (bc.second)),
                                                            std::get<2> (bc.second),
                                                            bc.first);
            }
        }
        
        //copy dirichlet boundary conditions
        dolfin::log (dolfin::DBG, "Copying Dirichlet boundary conditions...");
        for (auto& bc : this->dirichletBCs_)
        {
            clonedProblem->addDirichletBC (bc.second, bc.first);
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
    void TimeDependentProblem::advanceTime ()
    {
        t_ += dt_;
        
        dolfin::log (dolfin::INFO, "TIME = %f s", t_);

        std::string previousSolutionName = parameters ["previous_solution_name"];
        std::vector<std::string> previousSolutionCoefficientTypes;
        parameters ("previous_solution_coefficient_types").get_parameter_keys (previousSolutionCoefficientTypes);
        int timeSteppingSolutionComponent = parameters ["time_stepping_solution_component"];
        
        bool previousSolutionIsSetExternally = parameters ["previous_solution_is_set_externally"];
        if (previousSolutionIsSetExternally == false)
        {
            dolfin::begin (dolfin::DBG, "Setting previous solution coefficients...");
            for (unsigned int step = 1; step <= nTimeSchemeSteps_; ++step)
            {
                // if step > 1, append step number to the name of the coefficient to set 
                if (step > 1)
                {
                    previousSolutionName = previousSolutionName + "_" + std::to_string (step);
                }
                
                // get the index of the function in solution_ to use to set the coefficient
                unsigned int index = solution_.size () - step;
                for (auto& previousSolutionCoefficientType : previousSolutionCoefficientTypes)
                {
                    if (timeSteppingSolutionComponent >= 0)
                    {
                        timeSteppingProblem_ -> setCoefficient 
                            (previousSolutionCoefficientType, 
                             dolfin::reference_to_no_delete_pointer ((solution_ [index].second) [timeSteppingSolutionComponent]), 
                             previousSolutionName);
                    }
                    else
                    {
                        timeSteppingProblem_ -> setCoefficient 
                            (previousSolutionCoefficientType, 
                             dolfin::reference_to_no_delete_pointer (solution_ [index].second),
                             previousSolutionName);
                    }
                } 
            }
            dolfin::end ();
        }
        else
        {
            dolfin::log (dolfin::DBG, "Skipping previous solution setting loop since it is set externally.");
        }
    }
            
    

    void TimeDependentProblem::setTimeDependentDirichletBCs ()
    {
        dolfin::begin (dolfin::DBG, "Setting time dependent Dirichlet's boudary conditions...");
        
        // reset time dependent Dirichlet BCs
        for (auto bcIterator = timeDependentDirichletBCs_.begin (); 
             bcIterator != timeDependentDirichletBCs_.end (); 
             bcIterator++)
        {
            resetTimeDependentDirichletBC (bcIterator);
        }

        dolfin::end ();
    }

    
    
    void TimeDependentProblem::resetTimeDependentDirichletBC 
        (std::map <TimeDependentProblem::TimeDependentDirichletBCKey, 
                   TimeDependentProblem::TimeDependentDirichletBCValue>
                   ::iterator bcIterator)
    {
        // get bc name
        std::string bcName = bcIterator->first;
        
        // get dcp::TimeDependentExpression object and set time equal to t_
        std::shared_ptr<dcp::TimeDependentExpression> condition = std::get<0> (bcIterator->second);
        condition->setTime (t_);
        
        // get dcp::Subdomain object
        std::shared_ptr<dcp::Subdomain> boundary = std::get<1> (bcIterator->second);

        // get component on which to enforce the bc
        int component = std::get<2> (bcIterator->second);
        
        dolfin::begin (dolfin::DBG, 
                       "Resetting time dependent Dirichlet's boundary condition \"%s\"...", 
                       bcName.c_str ());
        
        // remove the bc from timeSteppingProblem_
        removeDirichletBC (bcName);
        
        // remember that component = -1 means that the bc should be enforced on all the function space components
        if (component == -1)
        {
            addDirichletBC (condition, boundary, bcName);
        }
        else
        {
            addDirichletBC (condition, boundary, component, bcName);
        } 
        
        dolfin::end ();
    }

    

    void TimeDependentProblem::setTimeDependentCoefficients ()
    {
        dolfin::begin (dolfin::DBG, "Setting time dependent coefficients...");
        
        for (auto& coefficientPair : timeDependentCoefficients_)
        {
            std::shared_ptr <dcp::TimeDependentExpression> expression (coefficientPair.second);
            std::string coefficientName = std::get<0> (coefficientPair.first);
            std::string coefficientType = std::get<1> (coefficientPair.first);
            dolfin::begin (dolfin::DBG, 
                         "Coefficient: name \"%s\", type \"%s\"", 
                         coefficientName.c_str (),
                         coefficientType.c_str ());
            
            dolfin::log (dolfin::DBG, "Setting time in time dependent expression...");
            expression->setTime (t_);
            
            timeSteppingProblem_->setCoefficient (coefficientType, expression, coefficientName);
            
            dolfin::end ();
        }
        
        dolfin::end ();
    }
    
            
    
    void TimeDependentProblem::printFinishedWarning ()
    {
        dolfin::warning ("No time iteration performed in solve() function. End time already reached.");
    }
    


    void TimeDependentProblem::storeSolution (const dolfin::Function& solution, 
                                              const int& timeStep, 
                                              const int& storeInterval)
    {
        if (storeInterval > 0 && timeStep % storeInterval == 0)
        {
            dolfin::log (dolfin::DBG, "Saving time stepping problem solution in solutions vector...");
            solution_.push_back (std::make_pair (t_, solution));
        }
    } 
    


    void TimeDependentProblem::storeLastStepSolution (const dolfin::Function& solution, 
                                                      const int& timeStep, 
                                                      const int& storeInterval)
    {
        if (!(storeInterval > 0 && timeStep % storeInterval == 0))
        {
            dolfin::log (dolfin::DBG, "Saving last time step solution in solutions vector...");
            solution_.push_back (std::make_pair (t_, solution));
        }
    }

    

    void TimeDependentProblem::plotSolution (dolfin::Function& solution, 
                                             const int& timeStep, 
                                             const int& plotInterval, 
                                             const int& plotComponent,
                                             const bool& pause)
    {
        if (plotInterval > 0 && timeStep % plotInterval == 0)
        {
            dolfin::begin (dolfin::DBG, "Plotting...");

            // auxiliary variable, to enhance readability
            std::shared_ptr<dolfin::Function> functionToPlot;
            
            std::string plotTitle = parameters ["plot_title"];

            // if plotTitle is not empty, we need to prepend ", " so that the plot title is readable
            if (!plotTitle.empty ())
            {
                plotTitle = ", " + plotTitle;
            }

            // get right function to plot
            if (plotComponent == -1)
            {
                functionToPlot = dolfin::reference_to_no_delete_pointer (solution);
                dolfin::log (dolfin::DBG, "Plotting time stepping problem solution, all components...");
            }
            else
            {
                functionToPlot = dolfin::reference_to_no_delete_pointer (solution [plotComponent]);
                dolfin::log (dolfin::DBG, "Plotting time stepping problem solution, component %d...", plotComponent);
            }

            // actual plotting
            if (solutionPlotter_ == nullptr)
            {
                dolfin::log (dolfin::DBG, "Plotting in new dolfin::VTKPlotter object...");
                solutionPlotter_ = dolfin::plot (functionToPlot, "Time = " + std::to_string (t_) + plotTitle);
            }
            else if (! solutionPlotter_ -> is_compatible (functionToPlot))
            {
                dolfin::log (dolfin::DBG, "Existing plotter is not compatible with object to be plotted.");
                dolfin::log (dolfin::DBG, "Creating new dolfin::VTKPlotter object...");
                solutionPlotter_ = dolfin::plot (functionToPlot, "Time = " + std::to_string (t_) + plotTitle);
            }
            else 
            {
                solutionPlotter_ -> parameters ["title"] = std::string ("Time = " + std::to_string (t_) + plotTitle);
                solutionPlotter_ -> plot (functionToPlot);
            }

            if (pause)
            {
                dolfin::interactive ();
            }

            dolfin::end ();
        }
    }
}
