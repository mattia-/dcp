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

#include <dcp/optimizers/TimeDependentBacktrackingOptimizer.h>
#include <dolfin/parameter/Parameters.h>
#include <dolfin/log/dolfin_log.h>
#include <dolfin/function/Function.h>
#include <dolfin/common/NoDeleter.h>
#include <dolfin/fem/assemble.h>
#include <dolfin/fem/Form.h>
#include <string>
#include <functional>
#include <cmath>
#include <iomanip>
#include <dcp/utils/dotproductforms.h>

namespace dcp
{
    TimeDependentBacktrackingOptimizer::TimeDependentBacktrackingOptimizer ():
        BacktrackingOptimizer ()
    {
        dolfin::begin (dolfin::DBG, "Creating TimeDependentBacktrackingOptimizer object...");

        dolfin::log (dolfin::DBG, "Setting up parameters...");
        parameters["descent_method"] = "time_dependent_backtracking_gradient_method";
        dolfin::log (dolfin::DBG, "TimeDependentBacktrackingOptimizer object created");

        dolfin::end ();
    }



    void TimeDependentBacktrackingOptimizer::apply 
        (const std::vector<std::shared_ptr<dcp::GenericEquationSystem> > systems,
         const dcp::GenericObjectiveFunctional& objectiveFunctional, 
         dolfin::Function& initialGuess,
         const dcp::GenericDescentMethod::Updater& updater)
    {
        // check if systems are dcp::TimeDependentEquationSystem
        std::shared_ptr<dcp::TimeDependentEquationSystem> 
            primalSystem = std::dynamic_pointer_cast<dcp::TimeDependentEquationSystem> (systems[0]);
        if (primalSystem == nullptr)
        {
            dolfin::dolfin_error 
                ("dcp: TimeDependentBacktrackingOptimizer.cpp",
                 "apply",
                 "Primal system (aka systems[0]) is not a dcp::TimeDependentEquationSystem");
        }

        std::shared_ptr<dcp::TimeDependentEquationSystem> 
            adjointSystem = std::dynamic_pointer_cast<dcp::TimeDependentEquationSystem> (systems[1]);
        if (adjointSystem == nullptr)
        {
            dolfin::dolfin_error 
                ("dcp: TimeDependentBacktrackingOptimizer.cpp",
                 "apply",
                 "Adjoint system (aka systems[1]) is not a dcp::TimeDependentEquationSystem");
        }

        // check if times are ok
        if (dolfin::near (primalSystem->startTime (), adjointSystem->endTime ()) == false)
        {
            dolfin::dolfin_error
                ("dcp: TimeDependentBacktrackingOptimizer.cpp",
                 "apply",
                 "Primal start time and adjoint end time mismatch");
        }
        
        if (dolfin::near (primalSystem->endTime (), adjointSystem->startTime ()) == false)
        {
            dolfin::dolfin_error
                ("dcp: TimeDependentBacktrackingOptimizer.cpp",
                 "apply",
                 "Primal end time and adjoint start time mismatch");
        }

        if (dolfin::near (primalSystem->dt (), -adjointSystem->dt ()) == false)
        {
            dolfin::dolfin_error
                ("dcp: TimeDependentBacktrackingOptimizer.cpp",
                 "apply",
                 "Primal dt and adjoint dt mismatch");
        }

        // set purge_inteval equal to 0, because we need all the solutions when we solve the backward-in-time adjoint
        // system 
        dolfin::begin (dolfin::DBG, "Setting parameter \"purge_interval\" to 0 for all problems...");
        for (std::size_t i = 0; i < systems.size (); ++i)
        {
            for (std::size_t j = 0; j < systems[i]->size (); ++j)
            {
                (*(systems[i]))[j].parameters["purge_interval"] = 0;
            }
        }
        dolfin::end ();
        
        // now revert back to base class apply method
        dcp::BacktrackingOptimizer::apply (systems, objectiveFunctional, initialGuess, updater);
    }


    // ************************************* //
    // ********** PRIVATE MEMBERS ********** //
    // ************************************* //
    void TimeDependentBacktrackingOptimizer::solve_ 
        (const std::vector<std::shared_ptr<dcp::GenericEquationSystem> > systems,
         const std::string& solveType)
    {
        // get primal and adjoint problems
        // use of static_cast is justified by the fact that apply method already checks if systems are of the right type
        dcp::TimeDependentEquationSystem& primalSystem = 
            static_cast<dcp::TimeDependentEquationSystem&> (*(systems[0]));
        dcp::TimeDependentEquationSystem& adjointSystem = 
            static_cast<dcp::TimeDependentEquationSystem&> (*(systems[1]));

        if (solveType != "all" && solveType != "primal" && solveType != "adjoint")
        {
            dolfin::dolfin_error ("dcp: TimeDependentBacktrackingOptimizer.cpp",
                                  "solve",
                                  "Unknown solve type \"%s\"", 
                                  solveType.c_str ());
        }

        // 1) solve primal problem
        if (solveType == "all" || solveType == "primal")
        {
            dolfin::begin (dolfin::PROGRESS, "Initializing primal system...");
            initializePrimalSystem_ (primalSystem);
            dolfin::end (); // Initializing primal system

            dolfin::begin (dolfin::PROGRESS, "Solving primal system...");
            solvePrimalSystem_ (primalSystem);
            dolfin::end (); // Solving primal system
        }

        // 2) solve adjoint problem
        if (solveType == "all" || solveType == "adjoint")
        {
            dolfin::begin (dolfin::PROGRESS, "Initializing adjoint system...");
            initializeAdjointSystem_ (adjointSystem, primalSystem);
            dolfin::end (); // Initializing adjoint system

            dolfin::begin (dolfin::PROGRESS, "Solving adjoint system...");
            solveAdjointSystem_ (adjointSystem, primalSystem);
            dolfin::end (); // Solving adjoint system
        }
    }



    void TimeDependentBacktrackingOptimizer::initializePrimalSystem_ (dcp::TimeDependentEquationSystem& primalSystem)
    {

        for (std::size_t i = 0; i < primalSystem.size (); ++i)
        {
            dolfin::begin (dolfin::DBG, "Initializing problem number %d in primal problem...", i);
            primalSystem[i].clear ();
            primalSystem[i].restoreState ("initial_state", true);
            dolfin::end (); // Initializing problem number %d in primal problem
        }
    }



    void TimeDependentBacktrackingOptimizer::initializeAdjointSystem_ 
        (dcp::TimeDependentEquationSystem& adjointSystem,
         const dcp::TimeDependentEquationSystem& primalSystem)
    {
        for (std::size_t i = 0; i < adjointSystem.size (); ++i)
        {
            dolfin::begin (dolfin::DBG, "Initializing problem number %d in adjoint problem...", i);
            adjointSystem[i].clear ();
            adjointSystem[i].restoreState ("initial_state", true);
            dolfin::end (); // Initializing problem number %d in adjoint problem
        }
    }



    void TimeDependentBacktrackingOptimizer::solvePrimalSystem_ (dcp::TimeDependentEquationSystem& primalSystem)
    {
        primalSystem.solve ();
    }



    void TimeDependentBacktrackingOptimizer::solveAdjointSystem_ (dcp::TimeDependentEquationSystem& adjointSystem,
                                                                  const dcp::TimeDependentEquationSystem& primalSystem)
    {
        // i) Reserve space for solutions vector of each problem
        for (std::size_t i = 0; i < adjointSystem.size (); ++i)
        {
            adjointSystem[i].reserve ();
        }

        // ii) Solutions loop
        dolfin::begin (dolfin::DBG, "Solution loop...");
        unsigned int timeStep = 0;
        while (adjointSystem.isFinished () == 0)
        {
            timeStep++;
            dolfin::begin (dolfin::PROGRESS, "===== Timestep %d =====", timeStep);

            adjointSolveStepPreprocessing_ (adjointSystem, primalSystem, timeStep);

            adjointSystem.solve ("step");

            // plot and write to file problems solutions
            for (const auto& problemName : adjointSystem.problemsNames ())
            {
                // plot solution
                dcp::TimeDependentProblem& problem = adjointSystem[problemName];
                int plotInterval = problem.parameters["plot_interval"];

                if (plotInterval > 0 && timeStep % plotInterval == 0)
                {
                    problem.plotSolution ("last");
                }

                // write solution to file 
                int writeInterval = problem.parameters ["write_interval"];
                if (writeInterval > 0 && timeStep % writeInterval == 0)
                {
                    problem.writeSolutionToFile ("last");
                }
            }

            dolfin::end (); // Timestep %d
        }

        dolfin::end (); // Solution loop
    }



    void TimeDependentBacktrackingOptimizer::primalSolveStepPreprocessing_
        (dcp::TimeDependentEquationSystem& primalSystem,
         const dcp::TimeDependentEquationSystem& adjointSystem,
         unsigned int timeStep)
    {
    }



    void TimeDependentBacktrackingOptimizer::adjointSolveStepPreprocessing_ 
        (dcp::TimeDependentEquationSystem& adjointSystem,
         const dcp::TimeDependentEquationSystem& primalSystem,
         unsigned int timeStep)
    {
    }



    void TimeDependentBacktrackingOptimizer::update_ 
        (const std::vector<std::shared_ptr<dcp::GenericEquationSystem> > systems,
         const dcp::GenericDescentMethod::Updater& updater,
         const dolfin::GenericFunction& control)
    {
        // updater is called on the first element of the systems vector, since it represents the primal system, and the
        // control can only be enforced on the primal system (otherwise it would not be the primal system, would it?)
        updater (*(systems[0]), control);
    }
}
