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

#include <dcp/problems/TimeDependentProblem.h>
#include <dolfin/log/dolfin_log.h>
#include <dolfin/io/File.h>
#include <regex>

namespace dcp
{
    /******************* CONSTRUCTORS *******************/
    TimeDependentProblem::TimeDependentProblem 
        (const std::shared_ptr<dcp::GenericProblem> timeSteppingProblem,
         const std::shared_ptr<dcp::Time> time,
         const double& startTime,
         const double& dt,
         const double& endTime,
         std::initializer_list<std::string> dtCoefficientTypes,
         std::initializer_list<std::string> previousSolutionCoefficientTypes,
         const unsigned int& nTimeSchemeSteps)
        : 
            GenericProblem (timeSteppingProblem->functionSpace ()),
            timeSteppingProblem_ (timeSteppingProblem),
            time_ (time),
            startTime_ (startTime),
            dt_ (dt),
            endTime_ (endTime),
            timeDependentCoefficients_ (),
            timeDependentDirichletBCs_ (),
            nTimeSchemeSteps_ (nTimeSchemeSteps),
            timeDependentDirichletBCsCounter_ (0)
    { 
        dolfin::begin (dolfin::DBG, "Building TimeDependentProblem...");
        
        time_ -> setTo (startTime_);
        
        // we need a solution for every step in the time scheme! 
        // And the time should be t0-(nTimeSchemeSteps-1)*dt, to-(nTimeSchemeSteps-2)*dtdt, ... , t0-dt (t0 is the start time)
        dolfin::begin (dolfin::DBG, "Creating initial solutions...");
        unsigned int stepNumber;
        for (stepNumber = nTimeSchemeSteps_; stepNumber > 0; --stepNumber)
        {
            solution_.emplace_back (std::make_pair (time_ -> value () - (stepNumber - 1) * dt_, 
                                                    dolfin::Function (timeSteppingProblem_->functionSpace ())));
        }
        dolfin::log (dolfin::DBG, "Created %d initial solutions", stepNumber + 1);
        
        dolfin::end ();
        
        dolfin::log (dolfin::DBG, "Setting up parameters...");
        parameters.add ("problem_type", "time_dependent");
        parameters.add ("dt_name", "dt");
        parameters.add ("previous_solution_name", "u_old");
        parameters.add ("write_interval", 0);
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



    TimeDependentProblem::TimeDependentProblem 
        (const std::shared_ptr<dcp::GenericProblem> timeSteppingProblem,
         const double& startTime,
         const double& dt,
         const double& endTime,
         std::initializer_list<std::string> dtCoefficientTypes,
         std::initializer_list<std::string> previousSolutionCoefficientTypes,
         const unsigned int& nTimeSchemeSteps)
        : 
            GenericProblem (timeSteppingProblem->functionSpace ()),
            timeSteppingProblem_ (timeSteppingProblem),
            time_ (new dcp::Time (startTime)),
            startTime_ (startTime),
            dt_ (dt),
            endTime_ (endTime),
            timeDependentCoefficients_ (),
            timeDependentDirichletBCs_ (),
            nTimeSchemeSteps_ (nTimeSchemeSteps),
            timeDependentDirichletBCsCounter_ (0)
    { 
        dolfin::begin (dolfin::DBG, "Building TimeDependentProblem...");
        
        // we need a solution for every step in the time scheme! 
        // And the time should be t0, t0+dt, t0+2dt ... (t0 is the start time)
        dolfin::begin (dolfin::DBG, "Creating initial solutions...");
        unsigned int stepNumber;
        for (stepNumber = 0; stepNumber < nTimeSchemeSteps_; ++stepNumber)
        {
            solution_.emplace_back (std::make_pair (time_ -> value () + stepNumber * dt, 
                                                    dolfin::Function (timeSteppingProblem_->functionSpace ())));
        }
        dolfin::log (dolfin::DBG, "Created %d initial solutions", stepNumber + 1);
        
        // set the correct value for time_, since it was not incremented during the previous loop.
        // Use stepNumber - 1 since time_ will be incremented at the *beginning* of the time loop, not at the end
        time_ -> add ((stepNumber - 1) * dt);
            
        dolfin::end ();
        
        dolfin::log (dolfin::DBG, "Setting up parameters...");
        parameters.add ("problem_type", "time_dependent");
        parameters.add ("dt_name", "dt");
        parameters.add ("previous_solution_name", "u_old");
        parameters.add ("write_interval", 0);
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
    


    const std::vector <std::pair <double, dolfin::Function> >& TimeDependentProblem::solutionsVector () const
    {
        return solution_;
    }
    


    std::shared_ptr<dcp::Time> TimeDependentProblem::time () const
    {
        return time_;
    }
    


    const double& TimeDependentProblem::startTime ()
    {
        return startTime_;
    }
    


    const double& TimeDependentProblem::dt ()
    {
        return dt_;
    }



    const double& TimeDependentProblem::endTime ()
    {
        return endTime_;
    }
    


    dcp::GenericProblem& TimeDependentProblem::timeSteppingProblem ()
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
        // We add it to the existing value of nErasedElements so that we know whether both removal operations went fine
        nErasedElements += removeDirichletBC (bcName);
        dolfin::end ();

        dolfin::end ();
        
        // remember: both removal were ok if nErasedElements was equal to 1 both times!
        return nErasedElements == 2? true : false;
    }



    bool TimeDependentProblem::addTimeDependentCoefficient (const std::string& coefficientType,
                                                            std::shared_ptr<dcp::TimeDependentExpression> expression,
                                                            const std::string& coefficientName)
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
        return (dt_ > 0) ? (time_ -> value () >= endTime_ - DOLFIN_EPS) : (time_ -> value () <= endTime_ + DOLFIN_EPS);
    }
    


    void TimeDependentProblem::clear ()
    {
        // reset time_
        time_ -> setTo (startTime_);
        
        // clear solutions vector
        solution_.clear ();
        solution_.emplace_back (std::make_pair (time_ -> value (), 
                                                dolfin::Function (timeSteppingProblem_->functionSpace ()))
                                );
        
        // reset time dependent Dirichlet BCs
        for (auto bcIterator = timeDependentDirichletBCs_.begin (); 
             bcIterator != timeDependentDirichletBCs_.end (); 
             bcIterator++)
        {
            resetTimeDependentDirichletBC (bcIterator);
        }
    }
    


    void TimeDependentProblem::advanceTime ()
    {
        time_ -> add (dt_);
    }



    void TimeDependentProblem::solve (const std::string& solveType) 
    {
        // check solve type
        if (solveType != "default" && 
            solveType != "step" && 
            solveType != "clear+default" && 
            solveType != "clear+step" &&
            solveType != "steady" &&
            solveType != "stash")
        {
            dolfin::dolfin_error ("dcp: TimeDependentProblem.cpp", 
                                  "solve",
                                  "Unknown solve type \"%s\" requested",
                                  solveType.c_str ());
        }
        
        dolfin::log (dolfin::DBG, "Selected solve type: %s", solveType.c_str ());
        
        // check if we are already at the end of time loop
        // check only when solveType is not "steady" or "stash", since these solve type do not increment time
        if (isFinished () && solveType != "steady" && solveType != "stash")
        {
            printFinishedWarning ();
            return;
        }
            
        // ---- Set dt value in time stepping problem ---- //
        dolfin::begin (dolfin::DBG, "Setting dt value in time dependent problem...");
        
        // get parameters' values
        dolfin::Constant dt (dt_);
        std::string dtName = parameters ["dt_name"];
        std::vector<std::string> dtCoefficientTypes;
        parameters ("dt_coefficient_types").get_parameter_keys (dtCoefficientTypes);
       
        for (auto& i : dtCoefficientTypes)
        {
            timeSteppingProblem_->setCoefficient (i, dolfin::reference_to_no_delete_pointer (dt), dtName);
        } 
        
        dolfin::end ();
        
        // call clear() if needed ()
        if (std::regex_match (solveType, std::regex (".*clear.*")))
        {
            clear ();
        }
        
        // call right method depending on solveType
        if (solveType == "default" || solveType == "clear+default")
        {
             solveLoop (); 
        }
        else if (solveType == "step" || solveType == "clear+step")
        {
             step ();
        }
        else if (solveType == "steady")
        {
             steadySolve ();
        }
        else if (solveType == "stash")
        {
             stashSolve ();
        }
    }
    
    

    void TimeDependentProblem::applyStashedSolution ()
    {
        // save stashed solution in solution_, only if the time of the currently last solution is not the same as the
        // current time value. Otherwise, delete the last solution first
        if (solution_.back ().first == time_ -> value ())
        {
            solution_.pop_back ();
        }
        solution_.push_back (std::make_pair (time_ -> value (), stashedSolution_));
        
        stashedSolution_ = dolfin::Function (*functionSpace_);
    }
            


    void TimeDependentProblem::plotSolution (const std::string& plotType)
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
        
        if (plotType == "default")
        {
            for (auto& timeSolutionPair : solution_)
            {
                double time = timeSolutionPair.first;
                
                // get right function to plot
                if (plotComponent == -1)
                {
                    functionToPlot = dolfin::reference_to_no_delete_pointer (timeSolutionPair.second);
                    dolfin::log (dolfin::DBG, "Plotting problem solution, all components...");
                }
                else
                {
                    functionToPlot = dolfin::reference_to_no_delete_pointer (timeSolutionPair.second [plotComponent]);
                    dolfin::log (dolfin::DBG, "Plotting problem solution, component %d...", plotComponent);
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
        }
        else if (plotType == "last")
        {
            // get last element in vector
            auto timeSolutionPair = solution_.back ();
            
            // get time
            double time = timeSolutionPair.first;
            
            // get right function to plot
            if (plotComponent == -1)
            {
                functionToPlot = dolfin::reference_to_no_delete_pointer (timeSolutionPair.second);
                dolfin::log (dolfin::DBG, "Plotting problem solution, all components...");
            }
            else
            {
                functionToPlot = dolfin::reference_to_no_delete_pointer (timeSolutionPair.second [plotComponent]);
                dolfin::log (dolfin::DBG, "Plotting problem solution, component %d...", plotComponent);
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
        else if (plotType == "stashed")
        {
            // get right function to plot
            if (plotComponent == -1)
            {
                functionToPlot = dolfin::reference_to_no_delete_pointer (stashedSolution_);
                dolfin::log (dolfin::DBG, "Plotting stashed problem solution, all components...");
            }
            else
            {
                functionToPlot = dolfin::reference_to_no_delete_pointer (stashedSolution_ [plotComponent]);
                dolfin::log (dolfin::DBG, "Plotting stashed solution, component %d...", plotComponent);
            }

            // actual plotting
            if (solutionPlotter_ == nullptr)
            {
                dolfin::log (dolfin::DBG, "Plotting in new dolfin::VTKPlotter object...");
                solutionPlotter_ = dolfin::plot (functionToPlot,
                                                 "Time = " + std::to_string (time_->value()) + plotTitle + " (stashed)");
            }
            else if (! solutionPlotter_ -> is_compatible (functionToPlot))
            {
                dolfin::log (dolfin::DBG, "Existing plotter is not compatible with object to be plotted.");
                dolfin::log (dolfin::DBG, "Creating new dolfin::VTKPlotter object...");
                solutionPlotter_ = dolfin::plot (functionToPlot,
                                                 "Time = " + std::to_string (time_->value()) + plotTitle + " (stashed)");
            }
            else 
            {
                solutionPlotter_ -> parameters ["title"] = 
                    std::string ("Time = " + std::to_string (time_ -> value ()) + plotTitle + " (stashed)");
                solutionPlotter_ -> plot (functionToPlot);
            }

            if (pause)
            {
                dolfin::interactive ();
            }

        }
        else
        {
            dolfin::warning ("Uknown plot type \"%s\". No plot performed", plotType.c_str ());
        }
            
        dolfin::end ();
    }
    


    void TimeDependentProblem::writeSolutionToFile (const std::string& writeType)
    {
        // check if writeType is known
        if (writeType != "default" && writeType != "last" && writeType != "stashed")
        {
            dolfin::warning ("Uknown write type \"%s\". No write performed", writeType.c_str ());
            return;
        }
        
        // check if solutionFileName_ and parameters["solution_file_name"] coincide. If not, change solutionWriter_
        if (solutionFileName_ != std::string (parameters["solution_file_name"]))
        {
            solutionWriter_.reset (new dolfin::File (parameters ["solution_file_name"]));
            solutionFileName_ = std::string (parameters["solution_file_name"]);
        }
        
        if (writeType == "default")
        {
            for (const auto& element : solution_)
            {
                (*solutionWriter_) << std::make_pair (&(element.second), element.first);
            }
        }
        else if (writeType == "last")
        {
            (*solutionWriter_) << std::make_pair (&(solution_.back ().second), solution_.back ().first);
        }
        else if (writeType == "stashed")
        {
            (*solutionWriter_) << std::make_pair (&(stashedSolution_), time_ -> value ());
        }
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
                new dcp::TimeDependentProblem (timeSteppingProblem_, time_, startTime_, dt_, endTime_, {}, {});
            clonedProblem->timeDependentDirichletBCs_ = this->timeDependentDirichletBCs_;
        }
        else if (cloneMethod == "deep_clone")
        {
            // note that we pass an empty initializer_list to the constructor as dtCoefficientTypes and 
            // previousSolutionCoefficientTypes, because they will be copied when the parameters are copied anyway
            clonedProblem = 
                new dcp::TimeDependentProblem (timeSteppingProblem_, startTime_, dt_, endTime_, {}, {});
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
        for (auto& bc : this->dirichletBCs_)
        {
            clonedProblem->addDirichletBC (bc.second, bc.first);
        }

        for (auto& bc : this->timeDependentDirichletBCs_)
        {
            clonedProblem->addTimeDependentDirichletBC (*(std::get<0> (bc.second)),
                                                        *(std::get<1> (bc.second)),
                                                        std::get<2> (bc.second),
                                                        bc.first);
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
    void TimeDependentProblem::steadySolve ()
    {
        setTimeDependentCoefficients ();

        setTimeDependentDirichletBCs ();

        setPreviousSolutionsCoefficients ();
        
        dolfin::begin (dolfin::PROGRESS, "Solving time stepping problem...");
        timeSteppingProblem_->solve ("default");
        dolfin::end ();
        
        // save new solution in solution_, only if the time of the currently last solution is not the same as the
        // current time value. Otherwise, delete the last solution first
        if (dolfin::near (solution_.back ().first, time_ -> value ()))
        {
            solution_.pop_back ();
        }
        solution_.push_back (std::make_pair (time_ -> value (), timeSteppingProblem_->solution ("default")));
        
        // TODO save only the number of soultions needed to advance in time (i.e. one if the method is a one-step time
        // scheme, two for 2-steps time scheme....)
        
        // set stashedSolution_ to be equal to the last computed solution, so that when solution() is called
        // from a subiterations loop it gets the right one. In the case of "default" solveType, indeed, the
        // solution returned should be the same for all solution types
        stashedSolution_ = solution_.back ().second;
    }

    void TimeDependentProblem::stashSolve ()
    {
        setTimeDependentCoefficients ();

        setTimeDependentDirichletBCs ();

        setPreviousSolutionsCoefficients ();
        
        dolfin::begin (dolfin::PROGRESS, "Solving time stepping problem...");
        timeSteppingProblem_->solve ("default");
        dolfin::end ();
        
        // save new solution in stashedSolution_
        stashedSolution_ = timeSteppingProblem_->solution ("default");
    }
    


    void TimeDependentProblem::step ()
    {
        advanceTime ();
        
        dolfin::log (dolfin::PROGRESS, "TIME = %f s", time_ -> value ());
        
        steadySolve ();
    }

    
    
    void TimeDependentProblem::solveLoop ()
    {
        // ---- Problem settings ---- //
        int writeInterval = parameters ["write_interval"];
        int plotInterval = parameters ["plot_interval"];
        int plotComponent = parameters ["plot_component"];
        bool pause = parameters ["pause"];


        // ---- Problem solution ---- //
        dolfin::begin (dolfin::DBG, "Solving time dependent problem...");
        
        
        // function used to plot and save to file the solution through the time loop
        dolfin::Function tmpSolution = solution_.back ().second;
        
        // start time loop
        int timeStep = 0;
        while (!isFinished ())
        {
            timeStep++;
            
            dolfin::begin (dolfin::PROGRESS, "===== Timestep %d =====", timeStep);
            
            step ();
            
            tmpSolution = timeSteppingProblem_->solution ();
            
            // save solution to file according to time step and write interval.
            if (writeInterval > 0 && timeStep % writeInterval == 0)
            {
                dolfin::log (dolfin::DBG, "Saving time stepping problem solution in solutions vector...");
                writeSolutionToFile ("last");
            }
            
            // plot solution according to time step and plot interval
            plotSolution (tmpSolution, timeStep, plotInterval, plotComponent, pause);
            
            dolfin::end ();
        }
        
        // At this point, we just need to make sure that the solution on the last iteration was saved even though
        // timeStep % writeInterval != 0 if writeInterval is greater than 0 (but it must not be saved twice!)
        if (writeInterval > 0 && timeStep % writeInterval != 0)
        // negation of the condition of the last if, so it is performed if the last if was not
        {
            dolfin::log (dolfin::DBG, "Saving last time step solution in solutions vector...");
            writeSolutionToFile ("last");
        }
        
        dolfin::end ();
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
        
        // get dcp::TimeDependentExpression object and set time equal to time_
        std::shared_ptr<const dcp::TimeDependentExpression> condition = std::get<0> (bcIterator->second);
        condition -> time () -> setTo (time_ -> value ());
        
        // get dcp::Subdomain object
        std::shared_ptr<const dcp::Subdomain> boundary = std::get<1> (bcIterator->second);

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
            std::shared_ptr <const dcp::TimeDependentExpression> expression (coefficientPair.second);
            std::string coefficientName = std::get<0> (coefficientPair.first);
            std::string coefficientType = std::get<1> (coefficientPair.first);
            dolfin::begin (dolfin::DBG, 
                         "Coefficient: name \"%s\", type \"%s\"", 
                         coefficientName.c_str (),
                         coefficientType.c_str ());
            
            dolfin::log (dolfin::DBG, "Setting time in time dependent expression...");
            expression -> time () -> setTo (time_ -> value ());
            
            timeSteppingProblem_ -> setCoefficient (coefficientType, expression, coefficientName);
            
            dolfin::end ();
        }
        
        dolfin::end ();
    }
    
    

    void TimeDependentProblem::setPreviousSolutionsCoefficients ()
    {
        std::vector<std::string> previousSolutionCoefficientTypes;
        parameters ("previous_solution_coefficient_types").get_parameter_keys (previousSolutionCoefficientTypes);
        int timeSteppingSolutionComponent = parameters ["time_stepping_solution_component"];
        
        bool previousSolutionIsSetExternally = parameters ["previous_solution_is_set_externally"];
        if (previousSolutionIsSetExternally == false)
        {
            dolfin::begin (dolfin::DBG, "Setting previous solution coefficients...");
            for (unsigned int step = 1; step <= nTimeSchemeSteps_; ++step)
            {
                std::string previousSolutionName = parameters ["previous_solution_name"];
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



    void TimeDependentProblem::printFinishedWarning ()
    {
        dolfin::warning ("No time iteration performed in solve() function. End time already reached.");
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
                solutionPlotter_ = dolfin::plot (functionToPlot, "Time = " + 
                                                                 std::to_string (time_ -> value ()) + 
                                                                 plotTitle);
            }
            else if (! solutionPlotter_ -> is_compatible (functionToPlot))
            {
                dolfin::log (dolfin::DBG, "Existing plotter is not compatible with object to be plotted.");
                dolfin::log (dolfin::DBG, "Creating new dolfin::VTKPlotter object...");
                solutionPlotter_ = dolfin::plot (functionToPlot, "Time = " + 
                                                                 std::to_string (time_ -> value ()) + 
                                                                 plotTitle);  
            }
            else 
            {
                solutionPlotter_ -> parameters ["title"] = std::string ("Time = " + 
                                                                        std::to_string (time_ -> value ()) 
                                                                        + plotTitle);
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
