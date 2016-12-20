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
         const unsigned int& nTimeSchemeSteps) :
            GenericProblem (timeSteppingProblem->functionSpace ()),
            timeSteppingProblem_ (timeSteppingProblem),
            time_ (time),
            startTime_ (startTime),
            dt_ (dt),
            endTime_ (endTime),
            timestep_ (0),
            timeDependentExpressionCoefficients_ (),
            timeDependentExpressionDirichletBCs_ (),
            timeDependentFunctionCoefficients_ (),
            timeDependentFunctionDirichletBCs_ (),
            nTimeSchemeSteps_ (nTimeSchemeSteps),
            timeDependentDirichletBCsCounter_ (0),
            states_ ()
    {
        dolfin::begin (dolfin::DBG, "Building TimeDependentProblem...");

        time_->setTo (startTime_);

        // we need a solution for every step in the time scheme!
        // And the time should be:
        //      t0-(nTimeSchemeSteps-1)*dt, to-(nTimeSchemeSteps-2)*dt, ... , t0
        // (t0 is the start time)
        dolfin::begin (dolfin::DBG, "Creating initial solutions...");
        for (int stepNumber = nTimeSchemeSteps_ - 1; stepNumber >= 0; --stepNumber)
        {
            solution_.emplace_back (std::make_pair (time_->value () - stepNumber * dt_,
                                                    dolfin::Function (timeSteppingProblem_->functionSpace ())));
        }
        dolfin::log (dolfin::DBG, "Created %d initial solutions", nTimeSchemeSteps_);

        dolfin::end (); // Creating initial solutions...

        dolfin::log (dolfin::DBG, "Setting up parameters...");
        parameters.add ("problem_type", "time_dependent");
        parameters.add ("dt_name", "dt");
        parameters.add ("previous_solution_name", "u_old");
        parameters.add ("write_interval", 0);
        parameters.add ("plot_interval", 0);
        parameters.add ("purge_interval", int (nTimeSchemeSteps_));
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
         const unsigned int& nTimeSchemeSteps) :
            TimeDependentProblem (timeSteppingProblem,
                                  std::make_shared<dcp::Time> (startTime),
                                  startTime,
                                  dt,
                                  endTime,
                                  dtCoefficientTypes,
                                  previousSolutionCoefficientTypes,
                                  nTimeSchemeSteps)
    {

    }



    /******************* GETTERS *******************/
    std::shared_ptr<const dolfin::Mesh> TimeDependentProblem::mesh () const
    {
        return timeSteppingProblem_->mesh ();
    }



    std::shared_ptr<const dolfin::FunctionSpace> TimeDependentProblem::functionSpace () const
    {
        return timeSteppingProblem_->functionSpace ();
    }



    const dolfin::DirichletBC& TimeDependentProblem::dirichletBC (const std::string& bcName) const
    {
        return timeSteppingProblem_->dirichletBC (bcName);
    }



    const std::map<std::string, dolfin::DirichletBC>& TimeDependentProblem::dirichletBCs () const
    {
        return timeSteppingProblem_->dirichletBCs ();
    }



    const dcp::TimeDependentFunction& TimeDependentProblem::solutionsVector () const
    {
        return solution_;
    }



    std::shared_ptr<const dcp::Time> TimeDependentProblem::time () const
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
            dolfin::warning ("initial solution not set. Requested time step number %d is greater than problem's scheme's number of time steps",
                             stepNumber);
            return;
        }

        auto solutionsIterator = solution_.rbegin ();

        // the solution we want to set is identified by stepNumber. If it is equal to 1, we must set the last element
        // in solution_, if it is equal to 2 the last but one element and so on. That's why we subtract 1 from
        // stepNumber in the next instruction
        (solutionsIterator + (stepNumber - 1))->second = initialSolution;
    }



    void TimeDependentProblem::setInitialSolution (const dolfin::Expression& initialSolution,
                                                   const unsigned int& stepNumber)
    {
        if (stepNumber > nTimeSchemeSteps_)
        {
            dolfin::warning ("initial solution not set. Requested time step number %d is greater than problem's scheme's number of time steps",
                             stepNumber);
            return;
        }

        auto solutionsIterator = solution_.rbegin ();

        // the solution we want to set is identified by stepNumber. If it is equal to 1, we must set the last element
        // in solution_, if it is equal to 2 the last but one element and so on. That's why we subtract 1 from
        // stepNumber in the next instruction
        (solutionsIterator + (stepNumber - 1))->second = initialSolution;
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

        // check if bc with the same name is present in the map of time-dependent functions and stop here in this
        // case
        if (timeDependentFunctionDirichletBCs_.find (bcName) != timeDependentFunctionDirichletBCs_.end ())
        {
            dolfin::warning
                ("DirichletBC object not inserted because key \"%s\" already in map of time-dependent functions",
                 bcName.c_str ());
            return false;
        }

        dolfin::begin
            (dolfin::DBG,
             "Adding dirichlet boundary condition to time dependent boundary condition expressions map with name \"%s\"...",
             bcName.c_str ());

        auto result = timeDependentExpressionDirichletBCs_.emplace
            (bcName,
             std::make_tuple (std::shared_ptr<dcp::TimeDependentExpression> (condition.clone ()),
                              std::shared_ptr<dcp::Subdomain> (boundary.clone ()),
                              -1)
             );

        bool bcWasAdded = result.second;

        if (bcWasAdded == false)
        {
            dolfin::warning ("DirichletBC object \"%s\" not inserted in map",
                             bcName.c_str ());
        }
        else
        {
            dolfin::begin
                (dolfin::DBG,
                 "Adding dirichlet boundary condition to time stepping problem boundary conditions map with name \"%s\"...",
                 bcName.c_str ());
            bcWasAdded = bcWasAdded && addDirichletBC (std::get<0> ((result.first)->second),
                                                       std::get<1> ((result.first)->second),
                                                       bcName);
            dolfin::end (); // Adding dirichlet boundary condition to time stepping problem boundary conditions map
        }

        dolfin::end (); // Adding dirichlet boundary condition to time dependent boundary condition expressions map

        return result.second;
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

        // check if bc with the same name is present in the map of time-dependent functions and stop here in this
        // case
        if (timeDependentFunctionDirichletBCs_.find (bcName) != timeDependentFunctionDirichletBCs_.end ())
        {
            dolfin::warning
                ("DirichletBC object not inserted because key \"%s\" already in map of time-dependent functions",
                 bcName.c_str ());
            return false;
        }

        dolfin::begin
            (dolfin::DBG,
             "Adding dirichlet boundary condition to time dependent boundary condition expressions map with name \"%s\"...",
             bcName.c_str ());

        auto result = timeDependentExpressionDirichletBCs_.emplace
            (bcName,
             std::make_tuple (std::shared_ptr<dcp::TimeDependentExpression> (condition.clone ()),
                              std::shared_ptr<dcp::Subdomain> (boundary.clone ()),
                              component)
             );

        bool bcWasAdded = result.second;

        if (bcWasAdded == false)
        {
            dolfin::warning ("DirichletBC object \"%s\" not inserted in map",
                             bcName.c_str ());
        }
        else
        {
            dolfin::begin
                (dolfin::DBG,
                 "Adding dirichlet boundary condition to time stepping problem boundary conditions map with name \"%s\"...",
                 bcName.c_str ());
            bcWasAdded = bcWasAdded && addDirichletBC (std::get<0> ((result.first)->second),
                                                       std::get<1> ((result.first)->second),
                                                       component,
                                                       bcName);
            dolfin::end (); // Adding dirichlet boundary condition to time stepping problem boundary conditions map
        }

        dolfin::end (); // Adding dirichlet boundary condition to time dependent boundary condition expressions map

        return result.second;
    }



    bool TimeDependentProblem::addTimeDependentDirichletBC
        (const std::shared_ptr<const dcp::TimeDependentExpression> condition,
         const std::shared_ptr<const dcp::Subdomain> boundary,
         std::string bcName)
    {
        return addTimeDependentDirichletBC (*condition, *boundary, bcName);
    }



    bool TimeDependentProblem::addTimeDependentDirichletBC
        (const std::shared_ptr<const dcp::TimeDependentExpression> condition,
         const std::shared_ptr<const dcp::Subdomain> boundary,
         const std::size_t& component,
         std::string bcName)
    {
        return addTimeDependentDirichletBC (*condition, *boundary, component, bcName);
    }



    bool TimeDependentProblem::addTimeDependentDirichletBC (const dcp::TimeDependentFunction& condition,
                                                            const dcp::Subdomain& boundary,
                                                            std::string bcName)
    {
        if (bcName.empty ())
        {
            bcName = "time_dependent_dirichlet_condition_" + std::to_string (timeDependentDirichletBCsCounter_);
            timeDependentDirichletBCsCounter_++;
        }

        // check if bc with the same name is present in the map of time-dependent expressions and stop here in this
        // case
        if (timeDependentExpressionDirichletBCs_.find (bcName) != timeDependentExpressionDirichletBCs_.end ())
        {
            dolfin::warning
                ("DirichletBC object not inserted because key \"%s\" already in map of time-dependent expression",
                 bcName.c_str ());
            return false;
        }

        dolfin::begin
            (dolfin::DBG,
             "Adding dirichlet boundary condition to time dependent boundary condition functions map with name \"%s\"...",
             bcName.c_str ());

        auto result = timeDependentFunctionDirichletBCs_.emplace
            (bcName,
             std::make_tuple (dolfin::reference_to_no_delete_pointer (condition),
                              std::shared_ptr<dcp::Subdomain> (boundary.clone ()),
                              -1)
             );

        bool bcWasAdded = result.second;

        if (bcWasAdded == false)
        {
            dolfin::warning ("DirichletBC object \"%s\" not inserted in map",
                             bcName.c_str ());
        }
        else
        {
            dolfin::begin
                (dolfin::DBG,
                 "Adding dirichlet boundary condition to time stepping problem boundary conditions map with name \"%s\"...",
                 bcName.c_str ());

            // get one element in condition, does not matter which one. The boundary condition will be reset when the
            // solve functions are called. This is just to add a condition with the given name to the underlying
            // time stepping problem
            std::size_t position = 0;
            const dcp::TimeDependentFunction::value_type& conditionElement =
                (*(std::get<0> ((result.first)->second)))[position];

            dolfin::log (dolfin::DBG,
                         "Time in time-dependent BC at position %d is %f",
                         position,
                         conditionElement.first);

            bcWasAdded =
                bcWasAdded && addDirichletBC (dolfin::reference_to_no_delete_pointer (conditionElement.second),
                                              std::get<1> ((result.first)->second),
                                              bcName);
            dolfin::end (); // Adding dirichlet boundary condition to time stepping problem boundary conditions map with name \"%s\"
        }

        dolfin::end (); // Adding dirichlet boundary condition to time dependent boundary condition functions map

        return result.second;
    }


    bool TimeDependentProblem::addTimeDependentDirichletBC (const dcp::TimeDependentFunction& condition,
                                                            const dcp::Subdomain& boundary,
                                                            const std::size_t& component,
                                                            std::string bcName)
    {
        if (bcName.empty ())
        {
            bcName = "time_dependent_dirichlet_condition_" + std::to_string (timeDependentDirichletBCsCounter_);
            timeDependentDirichletBCsCounter_++;
        }

        // check if bc with the same name is present in the map of time-dependent expressions and stop here in this
        // case
        if (timeDependentExpressionDirichletBCs_.find (bcName) != timeDependentExpressionDirichletBCs_.end ())
        {
            dolfin::warning
                ("DirichletBC object not inserted because key \"%s\" already in map of time-dependent expression",
                 bcName.c_str ());
            return false;
        }

        dolfin::begin
            (dolfin::DBG,
             "Adding dirichlet boundary condition to time dependent boundary condition functions map with name \"%s\"...",
             bcName.c_str ());

        auto result = timeDependentFunctionDirichletBCs_.emplace
            (bcName,
             std::make_tuple (dolfin::reference_to_no_delete_pointer (condition),
                              std::shared_ptr<dcp::Subdomain> (boundary.clone ()),
                              component)
             );

        bool bcWasAdded = result.second;

        if (bcWasAdded == false)
        {
            dolfin::warning ("DirichletBC object \"%s\" not inserted in map",
                             bcName.c_str ());
        }
        else
        {
            dolfin::begin
                (dolfin::DBG,
                 "Adding dirichlet boundary condition to time stepping problem boundary conditions map with name \"%s\"...",
                 bcName.c_str ());

            // get one element in condition, does not matter which one. The boundary condition will be reset when the
            // solve functions are called. This is just to add a condition with the given name to the underlying
            // time stepping problem
            std::size_t position = 0;
            const dcp::TimeDependentFunction::value_type& conditionElement =
                (*(std::get<0> ((result.first)->second)))[position];

            dolfin::log (dolfin::DBG,
                         "Time in time-dependent BC at position %d is %f",
                         position,
                         conditionElement.first);

            bcWasAdded =
                bcWasAdded && addDirichletBC (dolfin::reference_to_no_delete_pointer (conditionElement.second),
                                              std::get<1> ((result.first)->second),
                                              component,
                                              bcName);
            dolfin::end (); // Adding dirichlet boundary condition to time stepping problem boundary conditions map with name \"%s\"
        }

        dolfin::end (); // Adding dirichlet boundary condition to time dependent boundary condition functions map

        return result.second;
    }



    bool TimeDependentProblem::addTimeDependentDirichletBC
        (const std::shared_ptr<const dcp::TimeDependentFunction> condition,
         const std::shared_ptr<const dcp::Subdomain> boundary,
         std::string bcName)
    {
        return addTimeDependentDirichletBC (*condition, *boundary, bcName);
    }



    bool TimeDependentProblem::addTimeDependentDirichletBC
        (const std::shared_ptr<const dcp::TimeDependentFunction> condition,
         const std::shared_ptr<const dcp::Subdomain> boundary,
         const std::size_t& component,
         std::string bcName)
    {
        return addTimeDependentDirichletBC (*condition, *boundary, component, bcName);
    }



    bool TimeDependentProblem::removeTimeDependentDirichletBC (const std::string& bcName)
    {
        dolfin::begin (dolfin::DBG,
                       "Removing dirichlet boundary condition \"%s\" from time dependent boundary conditions map...",
                       bcName.c_str ());
        std::size_t nErasedElements = timeDependentExpressionDirichletBCs_.erase (bcName);
        nErasedElements += timeDependentFunctionDirichletBCs_.erase (bcName);

        if (nErasedElements == 0)
        {
            dolfin::warning ("Cannot remove dirichlet boundary condition \"%s\" because it was not found in any time-dependent BCs map",
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



    bool TimeDependentProblem::addTimeDependentCoefficient
        (const std::string& coefficientType,
         const std::shared_ptr<dcp::TimeDependentExpression> expression,
         const std::string& coefficientName)
    {
        auto coefficientID = std::make_pair (coefficientName, coefficientType);

        // check if coefficient with the same name exists in the map of time-dependent functions and stop here in this
        // case
        if (timeDependentFunctionCoefficients_.find (coefficientID) != timeDependentFunctionCoefficients_.end ())
        {
            dolfin::warning
                ("Time dependent coefficient expression not inserted because key (\"%s\", \"%s\") is already present in map of time-dependent functions",
                 coefficientName.c_str (),
                 coefficientType.c_str ());
            return false;
        }

        dolfin::begin (dolfin::DBG,
                       "Inserting time dependent coefficient expression in map with name \"%s\" and type \"%s\"...",
                       coefficientName.c_str (),
                       coefficientType.c_str ());

        auto result = timeDependentExpressionCoefficients_.insert (std::make_pair (coefficientID, expression));

        if (result.second == false)
        {
            dolfin::warning ("Time dependent coefficient with key (\"%s\", \"%s\") not inserted in map",
                             coefficientName.c_str (),
                             coefficientType.c_str ());
        }

        dolfin::end ();

        return result.second;
    }



    bool TimeDependentProblem::addTimeDependentCoefficient
        (const std::string& coefficientType,
         const std::shared_ptr<const dcp::TimeDependentFunction> function,
         const std::string& coefficientName)
    {
        auto coefficientID = std::make_pair (coefficientName, coefficientType);

        // check if coefficient with the same name exists in the map of time-dependent expressions and stop here in this
        // case
        if (timeDependentExpressionCoefficients_.find (coefficientID) != timeDependentExpressionCoefficients_.end ())
        {
            dolfin::warning
                ("Time dependent coefficient expression not inserted because key (\"%s\", \"%s\") is already present in map of time-dependent functions",
                 coefficientName.c_str (),
                 coefficientType.c_str ());
            return false;
        }

        dolfin::begin (dolfin::DBG,
                       "Inserting time dependent coefficient function in map with name \"%s\" and type \"%s\"...",
                       coefficientName.c_str (),
                       coefficientType.c_str ());

        auto result = timeDependentFunctionCoefficients_.insert (std::make_pair (coefficientID, function));

        if (result.second == false)
        {
            dolfin::warning ("Time dependent coefficient with key (\"%s\", \"%s\") not inserted in map",
                             coefficientName.c_str (),
                             coefficientType.c_str ());
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

        std::size_t nErasedElements = timeDependentExpressionCoefficients_.erase (coefficientID);
        nErasedElements += timeDependentFunctionCoefficients_.erase (coefficientID);

        if (nErasedElements == 0)
        {
            dolfin::warning ("Time dependent coefficient with key (\"%s\", \"%s\") not removed: key was not found in any time-dependent coefficients map",
                             coefficientName.c_str (),
                             coefficientType.c_str ());
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
        return (dt_ > 0) ?  (time_->value () >= endTime_) : (time_->value () <= endTime_);
    }



    void TimeDependentProblem::clear ()
    {
        // reset time_
        time_->setTo (startTime_);
        timestep_ = 0;

        // clear and reinitialize solutions vector
        solution_.clear ();

        dolfin::begin (dolfin::DBG, "Creating initial solutions...");
        for (int stepNumber = nTimeSchemeSteps_ - 1; stepNumber >= 0; --stepNumber)
        {
            solution_.emplace_back (std::make_pair (time_->value () - stepNumber * dt_,
                                                    dolfin::Function (timeSteppingProblem_->functionSpace ())));
        }
        dolfin::log (dolfin::DBG, "Created %d initial solutions", nTimeSchemeSteps_);

        dolfin::end (); // Creating initial solutions...

        // reset time dependent Dirichlet BCs (expressions)
        for (auto bcIterator = timeDependentExpressionDirichletBCs_.begin ();
             bcIterator != timeDependentExpressionDirichletBCs_.end ();
             bcIterator++)
        {
            resetTimeDependentDirichletBC_ (bcIterator);
        }

        // reset time dependent Dirichlet BCs (functions)
        for (auto bcIterator = timeDependentFunctionDirichletBCs_.begin ();
             bcIterator != timeDependentFunctionDirichletBCs_.end ();
             bcIterator++)
        {
            resetTimeDependentDirichletBC_ (bcIterator);
        }
    }



    bool TimeDependentProblem::saveState (const std::string& stateName)
    {
        dolfin::begin (dolfin::DBG, "Saving state with name \"%s\"...", stateName.c_str ());

        // create state
        dcp::TimeDependentProblem::TimeDependentProblemState newState;

        // save time value in state
        std::get<0> (newState) = time_->value();

        // save timestep in state
        std::get<1> (newState) = timestep_;

        // resize state's vector to contain needed number of functions
        std::get<2> (newState).resize (nTimeSchemeSteps_, dolfin::Function (timeSteppingProblem_->functionSpace ()));

        // copy solutions to state's vector
        for (int stepNumber = nTimeSchemeSteps_ - 1; stepNumber >= 0; --stepNumber)
        {
            std::get<2> (newState)[nTimeSchemeSteps_ - 1 - stepNumber] = solution_[solution_.size() - 1 - stepNumber].second;
        }

        // insert state in map and return success/fail result
        auto result = states_.insert (std::make_pair (stateName, newState));

        if (result.second == false)
        {
            dolfin::warning ("State \"%s\" not inserted in map. Maybe a state with the same name already existed?",
                             stateName.c_str ());
        }
        else
        {
            dolfin::log (dolfin::DBG,
                         "State saved with name \"%s\": time value %f, timestep %d, size of solutions vector: %d",
                         stateName.c_str (),
                         std::get<0> (newState),
                         std::get<1> (newState),
                         std::get<2> (newState).size ());

        }


        dolfin::end (); // Saving state with name %s

        return result.second;
    }



    bool TimeDependentProblem::restoreState (const std::string& stateName, const bool& clear)
    {
        dolfin::begin (dolfin::DBG, "Restoring state with name \"%s\"...", stateName.c_str ());

        // get state
        auto stateIterator = states_.find (stateName);

        if (stateIterator == states_.end())
        {
            dolfin::warning ("State \"%s\" not restored. No state with given name found in stored states map",
                             stateName.c_str ());

            dolfin::end (); // Restoring state with name %s
            return false;
        }

        // clear solutions vector if requested by the user
        if (clear)
        {
            solution_.clear ();
        }

        const dcp::TimeDependentProblem::TimeDependentProblemState& state = stateIterator->second;

        // reset time
        dolfin::begin (dolfin::DBG, "Resetting time...");
        time_->setTo (std::get<0> (state));
        dolfin::log (dolfin::DBG, "Time is now %f", time_->value());
        dolfin::end (); // Resetting time

        // reset timestep
        dolfin::begin (dolfin::DBG, "Resetting timestep...");
        timestep_ = std::get<1> (state);
        dolfin::log (dolfin::DBG, "Timestep is now %d", timestep_);
        dolfin::end (); // Resetting timestep

        // copy solutions from state's vector to solutions_
        dolfin::begin (dolfin::DBG, "Setting solutions...");
        for (int stepNumber = nTimeSchemeSteps_ - 1; stepNumber >= 0; --stepNumber)
        {
            solution_.emplace_back (std::make_pair (time_->value () - stepNumber * dt_,
                                                    std::get<2> (state)[nTimeSchemeSteps_ - 1 - stepNumber]));
        }
        dolfin::log (dolfin::DBG, "Set %d solutions", std::get<2> (state).size ());
        dolfin::end (); // Setting solutions...

        dolfin::end (); // Restoring state with name %s

        return true;
    }



    bool TimeDependentProblem::deleteState (const std::string& stateName)
    {
        dolfin::begin (dolfin::DBG, "Deleting state with name \"%s\"...", stateName.c_str ());

        // get state
        auto stateIterator = states_.find (stateName);

        if (stateIterator == states_.end())
        {
            dolfin::warning ("State \"%s\" not deleted. No state with given name found in stored states map",
                             stateName.c_str ());

            dolfin::end (); // Deleting state with name %s
            return false;
        }

        states_.erase (stateIterator);

        dolfin::end (); // Deleting state with name %s

        return true;
    }



    void TimeDependentProblem::reserve (std::size_t nElements)
    {
        dolfin::begin (dolfin::DBG, "Setting solutions vector capacity...");
        solution_.reserve (nElements);
        dolfin::log (dolfin::DBG, "Solutions vector capacity is now %d", solution_.capacity());
        dolfin::end(); // Setting solutions vector capacity
    }



    void TimeDependentProblem::reserve ()
    {
        // reserve space for solution_, so that as few reallocations as possible happen.
        // The capacity needed is equal to the purge interval (if it is greater than the number of steps in the time
        // scheme) or the number of time steps to be performed (otherwise). Granted, this does not take into account the
        // cases where nSteps < nTimeSchemeSteps_, but they are marginal and not a big deal anyway, since the
        // sizes/capacities involved for solution_ would be very small (say O(10))
        int purgeInterval = parameters["purge_interval"];

        // Estimate number of steps to be performed
        int nSteps = int ((endTime_ - startTime_) / dt_);
        int neededCapacity = purgeInterval >= nTimeSchemeSteps_ ? purgeInterval + 1 : nSteps;

        this->reserve (neededCapacity);
    }



    void TimeDependentProblem::advanceTime ()
    {
        // increment timestep
        timestep_++;

        // set new time value
        time_ -> setTo (startTime_ + timestep_ * dt_);
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
            dolfin::warning ("No time iteration performed in solve() function. End time already reached. Time is %d",
                             time_->value ());
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

        dolfin::end (); // Setting dt value in time dependent problem

        // call clear() if needed ()
        if (std::regex_match (solveType, std::regex (".*clear.*")))
        {
            clear ();
        }

        // call right method depending on solveType
        if (solveType == "default" || solveType == "clear+default")
        {
            solveLoop_ ();
        }
        else if (solveType == "step" || solveType == "clear+step")
        {
            step_ ();
        }
        else if (solveType == "steady")
        {
            steadySolve_ ();
        }
        else if (solveType == "stash")
        {
            stashSolve_ ();
        }
    }



    void TimeDependentProblem::applyStashedSolution ()
    {
        // save stashed solution in solution_, only if the time of the currently last solution is not the same as the
        // current time value. Otherwise, delete the last solution first
        if (solution_.back ().first == time_->value ())
        {
            solution_.pop_back ();
        }
        solution_.push_back (std::make_pair (time_->value (), stashedSolution_));

        purgeSolutionsVector_ ();

        stashedSolution_ = dolfin::Function (functionSpace_);
    }



    void TimeDependentProblem::plotSolution (const std::string& plotType)
    {
        // check if plotType is known
        if (plotType != "default" && plotType != "last" && plotType != "stashed")
        {
            dolfin::warning ("Unknown plot type \"%s\". No plot performed", plotType.c_str ());
            return;
        }

        bool pause = parameters ["pause"];

        dolfin::begin (dolfin::DBG, "Plotting...");

        // get vector of plot components
        std::vector<int> plotComponents;
        std::stringstream plotComponentsStream ((std::string (parameters ["plot_components"])));

        // auxiliary variable to push the stream values into the vector
        int component;
        while (plotComponentsStream >> component)
        {
            plotComponents.push_back (component);
        }

        // check if solutionPlotters_ has right size. If not, clear it and reinitialize it with null pointers
        if (solutionPlotters_.size() != plotComponents.size())
        {
            solutionPlotters_.clear();
            solutionPlotters_.resize (plotComponents.size(), nullptr);
        }

        // auxiliary variable, to enhance readability
        std::shared_ptr<dolfin::Function> functionToPlot;

        if (plotType == "default")
        {
            for (auto& timeSolutionPair : solution_)
            {
                double time = timeSolutionPair.first;

                for (std::size_t i = 0; i < plotComponents.size (); ++i)
                {
                    int component = plotComponents[i];

                    // get right function to plot
                    if (component == -1)
                    {
                        functionToPlot = dolfin::reference_to_no_delete_pointer (timeSolutionPair.second);
                        dolfin::log (dolfin::DBG, "Plotting problem solution, all components...");
                    }
                    else
                    {
                        functionToPlot = dolfin::reference_to_no_delete_pointer (timeSolutionPair.second [component]);
                        dolfin::log (dolfin::DBG, "Plotting problem solution, component %d...", component);
                    }

                    std::string plotTitle = parameters ["plot_title"];
                    if (plotTitle.empty())
                    {
                        plotTitle = "Time = " + std::to_string (time);
                    }
                    else
                    {
                        plotTitle = "Time = " + std::to_string (time) + ", " + plotTitle;
                    }

                    if (component != -1)
                    {
                        plotTitle += ", component " + std::to_string (component);
                    }

                    // actual plotting
                    plot_ (solutionPlotters_[i], functionToPlot, plotTitle);
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

            for (std::size_t i = 0; i < plotComponents.size (); ++i)
            {
                int component = plotComponents[i];

                // get right function to plot
                if (component == -1)
                {
                    functionToPlot = dolfin::reference_to_no_delete_pointer (timeSolutionPair.second);
                    dolfin::log (dolfin::DBG, "Plotting problem solution, all components...");
                }
                else
                {
                    functionToPlot = dolfin::reference_to_no_delete_pointer (timeSolutionPair.second [component]);
                    dolfin::log (dolfin::DBG, "Plotting problem solution, component %d...", component);
                }

                std::string plotTitle = parameters ["plot_title"];
                if (plotTitle.empty())
                {
                    plotTitle = "Time = " + std::to_string (time);
                }
                else
                {
                    plotTitle = "Time = " + std::to_string (time) + ", " + plotTitle;
                }

                if (component != -1)
                {
                    plotTitle += ", component " + std::to_string (component);
                }

                // actual plotting
                plot_ (solutionPlotters_[i], functionToPlot, plotTitle);
            }

            if (pause)
            {
                dolfin::interactive ();
            }

        }
        else // aka plotType == "stashed", otherwise we would have exited on the first check
        {
            for (std::size_t i = 0; i < plotComponents.size (); ++i)
            {
                int component = plotComponents[i];

                // get right function to plot
                if (component == -1)
                {
                    functionToPlot = dolfin::reference_to_no_delete_pointer (stashedSolution_);
                    dolfin::log (dolfin::DBG, "Plotting stashed problem solution, all components...");
                }
                else
                {
                    functionToPlot = dolfin::reference_to_no_delete_pointer (stashedSolution_ [component]);
                    dolfin::log (dolfin::DBG, "Plotting stashed solution, component %d...", component);
                }

                std::string plotTitle = parameters ["plot_title"];
                if (plotTitle.empty())
                {
                    plotTitle = "Time = " + std::to_string (time_->value());
                }
                else
                {
                    plotTitle = "Time = " + std::to_string (time_->value()) + ", " + plotTitle;
                }

                if (component != -1)
                {
                    plotTitle += ", component " + std::to_string (component);
                }

                plotTitle += " (stashed)";

                // actual plotting
                plot_ (solutionPlotters_[i], functionToPlot, plotTitle);
            }

            if (pause)
            {
                dolfin::interactive ();
            }
        }

        dolfin::end ();
    }



    void TimeDependentProblem::writeSolutionToFile (const std::string& writeType)
    {
        dolfin::begin (dolfin::DBG, "Saving solution to file...");

        // check if writeType is known
        if (writeType != "default" && writeType != "last" && writeType != "stashed")
        {
            dolfin::warning ("Unknown write type \"%s\". No write performed", writeType.c_str ());
            return;
        }

        // get vector of write components
        std::vector<int> writeComponents;
        std::stringstream writeComponentsStream ((std::string (parameters ["write_components"])));

        // auxiliary variable to push the stream values into the vector
        int component;
        while (writeComponentsStream >> component)
        {
            writeComponents.push_back (component);
        }

        // get parameters values
        bool parameterBinaryWrite = bool(parameters["binary_write"]);
        std::string parameterSolutionFileName = std::string (parameters["solution_file_name"]);

        // check if solutionFileName_ and parameters["solution_file_name"] coincide, if solutionWriters_ has right
        // size, if writeComponents_ contains the same values as writeComponents and if binary writers should be used
        bool hasFileNameMismatch = (solutionFileName_ != parameterSolutionFileName);
        bool hasBinaryMismatch = (binaryWrite_ != parameterBinaryWrite);
        bool hasWritersSizeMismatch = parameterBinaryWrite == false ?
                                        solutionWriters_.size() != writeComponents.size() :
                                        binarySolutionWriters_.size () != writeComponents.size ();
        bool hasWriteComponentsMismatch = (writeComponents_ != writeComponents);
        if (hasFileNameMismatch || hasBinaryMismatch || hasWritersSizeMismatch || hasWriteComponentsMismatch)
        {
            solutionFileName_ = parameterSolutionFileName;
            writeComponents_ = writeComponents;
            binaryWrite_ = parameterBinaryWrite;

            if (binaryWrite_ == false)
            {
                solutionWriters_.clear ();
                solutionWriters_.resize (writeComponents_.size(), nullptr);
            }
            else
            {
                binarySolutionWriters_.clear ();
                binarySolutionWriters_.resize (writeComponents_.size(), nullptr);
            }
        }

        // at this point, the protected variables are fine and the writers just need to be created

        // auxiliary variable, to enhance readability
        std::shared_ptr<dolfin::Function> functionToWrite;

        if (writeType == "default")
        {
            for (auto& timeSolutionPair : solution_)
            {
                double time = timeSolutionPair.first;

                for (std::size_t i = 0; i < writeComponents.size (); ++i)
                {
                    int component = writeComponents[i];

                    // get right function to write
                    if (component == -1)
                    {
                        functionToWrite = dolfin::reference_to_no_delete_pointer (timeSolutionPair.second);
                    }
                    else
                    {
                        functionToWrite = dolfin::reference_to_no_delete_pointer (timeSolutionPair.second [component]);
                    }

                    // file name that keeps track also of the component to be written to file
                    std::string filenameWithComponent = solutionFileName_;
                    if (component != -1)
                    {
                        // regex matching the extension, aka all the characters after the last dot (included)
                        std::regex extensionRegex ("(\\.[^.]*$)");

                        // add the component number before the extension ($n is the n-th backreference of the match)
                        filenameWithComponent = std::regex_replace (filenameWithComponent,
                                                                    extensionRegex,
                                                                    "_component" + std::to_string (component) + "$1");
                    }

                    // actual writing
                    if (binaryWrite_ == false)
                    {
                        write_ (solutionWriters_[i], functionToWrite, filenameWithComponent, time);
                    }
                    else
                    {
                        write_ (binarySolutionWriters_[i], functionToWrite, filenameWithComponent, time);
                    }
                }
            }
        }
        else if (writeType == "last")
        {
            // get last element in vector
            auto timeSolutionPair = solution_.back ();

            // get time
            double time = timeSolutionPair.first;

            for (std::size_t i = 0; i < writeComponents.size (); ++i)
            {
                int component = writeComponents[i];

                // get right function to write
                if (component == -1)
                {
                    functionToWrite = dolfin::reference_to_no_delete_pointer (timeSolutionPair.second);
                }
                else
                {
                    functionToWrite = dolfin::reference_to_no_delete_pointer (timeSolutionPair.second [component]);
                }

                // file name that keeps track also of the component to be written to file
                std::string filenameWithComponent = solutionFileName_;
                if (component != -1)
                {
                    // regex matching the extension, aka all the characters after the last dot (included)
                    std::regex extensionRegex ("(\\.[^.]*$)");

                    // add the component number before the extension ($n is the n-th backreference of the match)
                    filenameWithComponent = std::regex_replace (filenameWithComponent,
                                                                extensionRegex,
                                                                "_component" + std::to_string (component) + "$1");
                }

                // actual writing
                if (binaryWrite_ == false)
                {
                    write_ (solutionWriters_[i], functionToWrite, filenameWithComponent, time);
                }
                else
                {
                    write_ (binarySolutionWriters_[i], functionToWrite, filenameWithComponent, time);
                }
            }
        }
        else // aka writeType == "stashed", otherwise we would have exited on the first check
        {
            double time = time_->value();
            for (std::size_t i = 0; i < writeComponents.size (); ++i)
            {
                int component = writeComponents[i];

                // get right function to write
                if (component == -1)
                {
                    functionToWrite = dolfin::reference_to_no_delete_pointer (stashedSolution_);
                }
                else
                {
                    functionToWrite = dolfin::reference_to_no_delete_pointer (stashedSolution_ [component]);
                }

                // file name that keeps track also of the component to be written to file
                std::string filenameWithComponent = solutionFileName_;
                if (component != -1)
                {
                    // regex matching the extension, aka all the characters after the last dot (included)
                    std::regex extensionRegex ("(\\.[^.]*$)");

                    // add the component number before the extension ($n is the n-th backreference of the match)
                    filenameWithComponent = std::regex_replace (filenameWithComponent,
                                                                extensionRegex,
                                                                "_component" + std::to_string (component) + "$1");
                }

                // actual writing
                if (binaryWrite_ == false)
                {
                    write_ (solutionWriters_[i], functionToWrite, filenameWithComponent, time);
                }
                else
                {
                    write_ (binarySolutionWriters_[i], functionToWrite, filenameWithComponent, time);
                }
            }
        }

        dolfin::end (); // Saving solution to file
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
            clonedProblem->timeDependentExpressionDirichletBCs_ = this->timeDependentExpressionDirichletBCs_;
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

        for (auto& bc : this->timeDependentExpressionDirichletBCs_)
        {
            clonedProblem->addTimeDependentDirichletBC (*(std::get<0> (bc.second)),
                                                        *(std::get<1> (bc.second)),
                                                        std::get<2> (bc.second),
                                                        bc.first);
        }

        for (auto& bc : this->timeDependentFunctionDirichletBCs_)
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
    void TimeDependentProblem::steadySolve_ ()
    {
        setTimeDependentCoefficients_ ();

        setTimeDependentDirichletBCs_ ();

        setPreviousSolutionsCoefficients_ ();

        dolfin::begin (dolfin::DBG, "Solving time stepping problem...");
        timeSteppingProblem_->solve ("default");
        dolfin::end (); // Solving time stepping problem...

        // save new solution in solution_, only if the time of the currently last solution is not the same as the
        // current time value. Otherwise, delete the last solution first
        if (dolfin::near (solution_.back ().first, time_->value ()))
        {
            solution_.pop_back ();
        }
        solution_.push_back (std::make_pair (time_->value (), timeSteppingProblem_->solution ("default")));

        purgeSolutionsVector_ ();

        // set stashedSolution_ to be equal to the last computed solution, so that when solution() is called
        // from a subiterations loop it gets the right one. In the case of "default" solveType, indeed, the
        // solution returned should be the same for all solution types
        stashedSolution_ = solution_.back ().second;
    }

    void TimeDependentProblem::stashSolve_ ()
    {
        setTimeDependentCoefficients_ ();

        setTimeDependentDirichletBCs_ ();

        setPreviousSolutionsCoefficients_ ();

        dolfin::begin (dolfin::DBG, "Solving time stepping problem...");
        timeSteppingProblem_->solve ("default");
        dolfin::end (); // Solving time stepping problem...

        // save new solution in stashedSolution_
        stashedSolution_ = timeSteppingProblem_->solution ("default");
    }



    void TimeDependentProblem::step_ ()
    {
        advanceTime ();

        dolfin::log (dolfin::PROGRESS, "TIME = %f s", time_->value ());

        steadySolve_ ();
    }



    void TimeDependentProblem::solveLoop_ ()
    {
        // reserve space for solutions vector
        this->reserve();

        // ---- Problem settings ---- //
        int writeInterval = parameters ["write_interval"];
        int plotInterval = parameters ["plot_interval"];

        // ---- Problem solution ---- //
        dolfin::begin (dolfin::PROGRESS, "Start solution loop...");

        // start time loop
        while (!isFinished ())
        {
            dolfin::begin (dolfin::PROGRESS, "===== Timestep %d =====", timestep_ + 1);

            step_ ();

            // save solution to file according to time step and write interval.
            if (writeInterval > 0 && timestep_ % writeInterval == 0)
            {
                writeSolutionToFile ("last");
            }

            // plot solution according to time step and plot interval
            if (plotInterval > 0 && timestep_ % plotInterval == 0)
            {
                plotSolution ("last");
            }

            dolfin::end (); // Timestep %d
        }

        // At this point, we just need to make sure that the solution on the last iteration was saved even though
        // timestep_ % writeInterval != 0 if writeInterval is greater than 0 (but it must not be saved twice!)
        if (writeInterval > 0 && timestep_ % writeInterval != 0)
        // negation of the condition of the last if, so it is performed if the last if was not
        {
            dolfin::log (dolfin::DBG, "Saving last time step solution in solutions vector...");
            writeSolutionToFile ("last");
        }

        dolfin::end (); // Start solution loop...
    }



    void TimeDependentProblem::setTimeDependentDirichletBCs_ ()
    {
        dolfin::begin (dolfin::DBG, "Setting time dependent Dirichlet's boundary conditions...");

        // reset time dependent Dirichlet BCs (expressions)
        for (auto bcIterator = timeDependentExpressionDirichletBCs_.begin ();
             bcIterator != timeDependentExpressionDirichletBCs_.end ();
             bcIterator++)
        {
            resetTimeDependentDirichletBC_ (bcIterator);
        }

        // reset time dependent Dirichlet BCs (functions)
        for (auto bcIterator = timeDependentFunctionDirichletBCs_.begin ();
             bcIterator != timeDependentFunctionDirichletBCs_.end ();
             bcIterator++)
        {
            resetTimeDependentDirichletBC_ (bcIterator);
        }

        dolfin::end (); // Setting time dependent Dirichlet's boundary conditions
    }



    void TimeDependentProblem::resetTimeDependentDirichletBC_
        (std::map <dcp::TimeDependentProblem::TimeDependentDirichletBCKey,
                   dcp::TimeDependentProblem::TimeDependentExpressionDirichletBCValue>
                   ::iterator bcIterator)
    {
        // get bc name
        std::string bcName = bcIterator->first;

        // get dcp::TimeDependentExpression object and set time equal to time_
        std::shared_ptr<dcp::TimeDependentExpression> condition = std::get<0> (bcIterator->second);
        condition->setTime (time_->value ());

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

        dolfin::end (); // Resetting time dependent Dirichlet's boundary condition %s
    }



    void TimeDependentProblem::resetTimeDependentDirichletBC_
        (std::map <dcp::TimeDependentProblem::TimeDependentDirichletBCKey,
                   dcp::TimeDependentProblem::TimeDependentFunctionDirichletBCValue>
                   ::iterator bcIterator)
    {
        // get bc name
        std::string bcName = bcIterator->first;

        // get right element in dcp::TimeDependentFunction object
        std::size_t position = dt_ > 0 ? timestep_ : (std::get<0> (bcIterator->second))->size () - timestep_;
        const dcp::TimeDependentFunction::value_type& condition = (*(std::get<0> (bcIterator->second)))[position];

        dolfin::log (dolfin::DBG,
                     "Time in time-dependent BC at position %d is %f",
                     position,
                     condition.first);

        // get dcp::Subdomain object
        std::shared_ptr<const dcp::Subdomain> boundary = std::get<1> (bcIterator->second);

        // get component on which to enforce the bc
        int component = std::get<2> (bcIterator->second);

        dolfin::begin (dolfin::DBG,
                       "Resetting time dependent Dirichlet's boundary condition \"%s\" (position in data vector: %d)...",
                       bcName.c_str (),
                       position);

        // remove the bc from timeSteppingProblem_
        removeDirichletBC (bcName);

        // remember that component = -1 means that the bc should be enforced on all the function space components
        if (component == -1)
        {
            addDirichletBC (dolfin::reference_to_no_delete_pointer (condition.second),
                            boundary,
                            bcName);
        }
        else
        {
            addDirichletBC (dolfin::reference_to_no_delete_pointer (condition.second),
                            boundary,
                            component,
                            bcName);
        }

        dolfin::end (); // Resetting time dependent Dirichlet's boundary condition %s
    }



    void TimeDependentProblem::setTimeDependentCoefficients_ ()
    {
        dolfin::begin (dolfin::DBG, "Setting time dependent coefficients (expressions)...");

        for (auto& coefficientPair : timeDependentExpressionCoefficients_)
        {
            std::shared_ptr <dcp::TimeDependentExpression> expression (coefficientPair.second);
            std::string coefficientName = coefficientPair.first.first;
            std::string coefficientType = coefficientPair.first.second;
            dolfin::begin (dolfin::DBG,
                           "Coefficient: name: \"%s\", type: \"%s\"",
                           coefficientName.c_str (),
                           coefficientType.c_str ());

            dolfin::log (dolfin::DBG, "Setting time in time dependent expression...");
            expression->setTime (time_->value ());

            timeSteppingProblem_->setCoefficient (coefficientType, expression, coefficientName);

            dolfin::end (); // Coefficient: name %s, type %s
        }
        dolfin::end (); // Setting time dependent coefficients (expressions)


        dolfin::begin (dolfin::DBG, "Setting time dependent coefficients (functions)...");
        for (auto& coefficientPair : timeDependentFunctionCoefficients_)
        {
            std::string coefficientName = coefficientPair.first.first;
            std::string coefficientType = coefficientPair.first.second;
            std::size_t position = dt_ > 0 ? timestep_ : coefficientPair.second->size () - timestep_;

            dolfin::begin (dolfin::DBG,
                           "Coefficient: name: \"%s\", type: \"%s\", position in data vector: %d",
                           coefficientName.c_str (),
                           coefficientType.c_str (),
                           position);

            auto& function = (*(coefficientPair.second))[position];
            dolfin::log (dolfin::DBG,
                         "Time in time-dependent coefficient at position %d is %f",
                         position,
                         function.first);

            timeSteppingProblem_->setCoefficient (coefficientType,
                                                  dolfin::reference_to_no_delete_pointer (function.second),
                                                  coefficientName);

            dolfin::end (); // Coefficient: name %s, type %s
        }

        dolfin::end (); // Setting time dependent coefficients (functions)
    }



    void TimeDependentProblem::setPreviousSolutionsCoefficients_ ()
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
                        timeSteppingProblem_->setCoefficient
                            (previousSolutionCoefficientType,
                             dolfin::reference_to_no_delete_pointer ((solution_ [index].second) [timeSteppingSolutionComponent]),
                             previousSolutionName);
                    }
                    else
                    {
                        timeSteppingProblem_->setCoefficient
                            (previousSolutionCoefficientType,
                             dolfin::reference_to_no_delete_pointer (solution_ [index].second),
                             previousSolutionName);
                    }
                }
            }
            dolfin::end (); // Setting previous solution coefficients...
        }
        else
        {
            dolfin::log (dolfin::DBG, "Skipping previous solution setting loop since it is set externally.");
        }
    }



    void TimeDependentProblem::purgeSolutionsVector_ ()
    {
        int purgeInterval = parameters["purge_interval"];

        while (purgeInterval >= nTimeSchemeSteps_ && purgeInterval < solution_.size ())
        {
            dolfin::log (dolfin::DBG, "Purging solution vector...");
            solution_.erase (solution_.begin ());
        }
    }



    void TimeDependentProblem::write_ (std::shared_ptr<dolfin::File>& writer,
                                       const std::shared_ptr<const dolfin::Function> function,
                                       const std::string& filename,
                                       const double& t)
    {
        if (writer == nullptr)
        {
            writer.reset (new dolfin::File (filename));
        }

        // try-catch block because we don't know if the file format allows std::pair on the operator <<
        try
        {
            (*writer) << std::make_pair (function.get(), t);
        }
        catch (std::runtime_error& e)
        {
            (*writer) << (*function);
        }
    }



    void TimeDependentProblem::write_ (std::shared_ptr<dolfin::HDF5File>& writer,
                                       const std::shared_ptr<const dolfin::Function> function,
                                       const std::string& filename,
                                       const double& t)
    {
        if (writer == nullptr)
        {
            writer.reset (new dolfin::HDF5File (MPI_COMM_WORLD, filename, "w"));
        }

        writer->write (*function, "solution", t);
    }
}
