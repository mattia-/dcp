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
    TimeDependentProblem::TimeDependentProblem (const std::shared_ptr<dcp::AbstractProblem> timeSteppingProblem,
                                                const double& startTime,
                                                const double& dt,
                                                const double& endTime,
                                                const std::vector<std::string>& dtCoefficientTypes,
                                                const std::vector<std::string>& previousSolutionCoefficientTypes,
                                                const int& storeInterval,
                                                const int& plotInterval,
                                                const std::string& dtName,
                                                const std::string& previousSolutionName) : 
            AbstractProblem (timeSteppingProblem->mesh (), timeSteppingProblem->functionSpace ()),
            timeSteppingProblem_ (timeSteppingProblem),
            t_ (startTime)
    { 
        dolfin::begin (dolfin::DBG, "Building TimeDependentProblem...");
        
        solution_.emplace_back (dolfin::Function (timeSteppingProblem->functionSpace ()));
            
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
        
        dolfin::log (dolfin::DBG, "TimeDependentProblem object created");
    }



    /******************* GETTERS *******************/
    std::shared_ptr<dolfin::Mesh> TimeDependentProblem::mesh () const
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



    bool TimeDependentProblem::addDirichletBC (const dolfin::DirichletBC& dirichletCondition, std::string bcName)
    {
        return timeSteppingProblem_->addDirichletBC (dirichletCondition, bcName);
    }



    bool TimeDependentProblem::addDirichletBC (dolfin::DirichletBC&& dirichletCondition, std::string bcName)
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
    bool TimeDependentProblem::isEnded ()
    {
        double endTime = parameters ["end_time"];
        return t_ >= endTime;
    }
    


    void TimeDependentProblem::clear ()
    {
        solution_.clear ();
        t_ = parameters ["start_time"];
    }
    


    void TimeDependentProblem::step ()
    {
        solve ("step");
    }
    

    
    void TimeDependentProblem::solve (const std::string& type) 
    {
        if (std::regex_match (type, std::regex (".*clear.*")))
        {
            clear ();
        }
        
        // ---- Problem settings ---- //
        dolfin::begin (dolfin::DBG, "Setting up time dependent problem...");
        
        bool stepFlag = std::regex_match (type, std::regex (".*step.*"));
        
        // get parameters' values
        double endTime = parameters ["end_time"];
        if (t_ >= endTime)
        {
            dolfin::warning ("No time iteration performed in solve() function. Protected paramter t >= end_time\n",
                             "t = %f\n",
                             "end_time = %f\n",
                             t_,
                             endTime);
            
            return;
        }
        
        dolfin::Constant dt (parameters ["dt"]);
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
        
        dolfin::end ();
        
        
        // ---- Problem solution ---- //
        dolfin::begin (dolfin::DBG, "Solving time dependent problem...");
        
        int timeStep = 0;
        
        // maxTimeSteps is used to check if the time loop should be ended at a certain timeStep.
        // It is set by default to -1, so that the loop never stops because it has reached the maximum number
        // of iteration. It is changed to 1, though, if stepFlag is true, that is if we only want to perform one
        // time step. In this case the time loop must (and will) exit after the first iteration
        int maxTimeSteps = -1;
        if (stepFlag)
        {
            maxTimeSteps = 1;
        }
        
        while (t_ < endTime + DOLFIN_EPS && timeStep != maxTimeSteps)
        {
            timeStep++;
            t_ += dt;
            
            dolfin::begin (dolfin::INFO, "===== Time = %f s =====", t_);
            dolfin::log (dolfin::DBG, "===== Timestep %d =====", timeStep);
            
            dolfin::log (dolfin::DBG, "Setting previous solution coefficients...");
            for (auto& i : previousSolutionCoefficientTypes)
            {
                if (timeSteppingSolutionComponent >= 0)
                {
                    timeSteppingProblem_->setCoefficient 
                        (i, 
                         dolfin::reference_to_no_delete_pointer (solution_.back () [timeSteppingSolutionComponent]), 
                         previousSolutionName);
                }
                else
                {
                    timeSteppingProblem_->setCoefficient 
                        (i, 
                         dolfin::reference_to_no_delete_pointer (solution_.back ()), 
                         previousSolutionName);
                }
            } 
            
            dolfin::begin (dolfin::INFO, "Solving time stepping problem...");
            timeSteppingProblem_->solve ();
            dolfin::end ();
            
            solution_.back () = timeSteppingProblem_->solution ();
            
            // save solution in solution_ according to time step and store interval.
            // Note that the current solution will always be the last element of solution_.
            // If the current solution is not to be stored, on the next iteration solution_.back () will be used as 
            // the previous iteration and then overwritten.
            // If the current solution is to be stored, we push_back it, so that we have two copies of it in the
            // vector and on the next iteration only the second one will be overwritten.
            if (storeInterval > 0 && timeStep % storeInterval == 0)
            {
                dolfin::log (dolfin::DBG, "Saving time stepping problem solution in solutions vector...");
                solution_.push_back (solution_.back ());
            }
            
            if (plotInterval > 0 && timeStep % plotInterval == 0)
            {
                dolfin::log (dolfin::DBG, "Plotting time stepping problem solution...");
                dolfin::plot (solution_.back ());
                
                if (pause)
                {
                    dolfin::interactive ();
                }
            }
            
            dolfin::end ();
        }
        
        // At this point, the solution on the last iteration will be already in solution_, automatically.
        // We just need to make sure that it isn't in there twice
        if (storeInterval > 0 && timeStep % storeInterval == 0)
        {
            dolfin::log (dolfin::DBG, "Removing time stepping problem solution (saved twice in solutions vector)...");
            solution_.pop_back ();
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
            std::vector<std::string> dtCoefficientTypes;
            parameters ("dt_coefficient_types").get_parameter_keys (dtCoefficientTypes);
            std::vector<std::string> previousSolutionCoefficientTypes;
            parameters ("previous_solution_coefficient_types").get_parameter_keys (previousSolutionCoefficientTypes);
            clonedProblem = 
                new dcp::TimeDependentProblem (this->timeSteppingProblem_,
                                               this->parameters ["start_time"],
                                               this->parameters ["dt"],
                                               this->parameters ["end_time"],
                                               dtCoefficientTypes,
                                               previousSolutionCoefficientTypes,
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
                new dcp::TimeDependentProblem (this->timeSteppingProblem_,
                                               this->parameters ["start_time"],
                                               this->parameters ["dt"],
                                               this->parameters ["end_time"],
                                               dtCoefficientTypes,
                                               previousSolutionCoefficientTypes,
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
        
        dolfin::end ();
        
        return clonedProblem;
    }
}
