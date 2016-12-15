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

#ifndef SRC_OPTIMIZERS_TIMEDEPENDENTBACKTRACKINGIMPLEMENTER_H_INCLUDE_GUARD
#define SRC_OPTIMIZERS_TIMEDEPENDENTBACKTRACKINGIMPLEMENTER_H_INCLUDE_GUARD

#include <dolfin/log/dolfin_log.h>
#include <dcp/problems/GenericProblem.h>
#include <dcp/problems/TimeDependentEquationSystem.h>
#include <dcp/optimizers/GenericImplementer.h>

namespace dcp
{
    /*! \class TimeDependentBacktrackingImplementer TimeDependentBacktrackingImplementer.h
     *  \brief Class that implements the specific methods needed by the backtracking algorithm in the time-dependent
     *  case. Derives from \c TimeDependentBacktrackingImplementer .
     *
     */
    template <class T_ControlVariable>
    class TimeDependentBacktrackingImplementer : public GenericImplementer<T_ControlVariable>
    {
        // ---------------------------------------------------------------------------------------------//

        public:
            /************************* CONSTRUCTORS ********************/
            //! Constructor [1]
            /*!
             *  \param updater the updater to be used to update the primal-adjoint system during the backtracking
             *  algorithm. Must be copy-constructible
             */
            TimeDependentBacktrackingImplementer
                    (const typename dcp::GenericImplementer<T_ControlVariable>::Updater& updater);

            //! Constructor [2]
            /*!
             *  \param updater the updater to be used to update the primal-adjoint system during the backtracking
             *  algorithm. Must be copy-constructible
             *  \param searchDirectionComputer the object used to compute the search direction. Must be
             *  copy-constructible
             */
            TimeDependentBacktrackingImplementer
                    (const typename dcp::GenericImplementer<T_ControlVariable>::Updater& updater,
                     const typename dcp::GenericImplementer<T_ControlVariable>::SearchDirectionComputer&
                            searchDirectionComputer);


            /************************* DESTRUCTOR ********************/
            //! Destructor
            virtual ~TimeDependentBacktrackingImplementer () {};


            /********************** METHODS ***********************/
            //! Solve the equation systems representing the primal and the adjoint problem.
            /*!
             *  This function calls the subfunctions \c initializePrimalSystem , \c initializeAdjointSystem ,
             *  \c primalSolveStepPreprocessing , \c adjointSolveStepPreprocessing , \c solvePrimalSystem and
             *  \c solveAdjointSystem . They all have a pretty standard and unrefined behaviour, since
             *  automatizing the treatment of all possible cases would be hard and lead to unintuitive API. If one needs
             *  a specific behaviour, they need to fine tune the class by deriving it and overriding the subfunctions
             *  called by this method.
             *
             *  NB: this function assumes that the systems passed in the arugment \c systems are of type
             *  \c dcp::TimeDependentEquationSystem
             *
             *
             *  \param systems the set of systems to be solved
             *  \param solveType the type of solve requested; possible values in this class:
             *  \li \c all
             *  \li \c primal
             *  \li \c adjoint
             *  with obvious meaning
             */
            virtual void solve (const std::vector<std::shared_ptr<dcp::GenericEquationSystem> > systems,
                                const std::string& solveType) override;

            //! Initialize primal problem for solving
            /*!
             *  This method relies on \c dcp::TimeDependentProblem::restoreState() method, and it assumes that the
             *  initial state to be restored is called \c "initial_state"
             *
             *  \param primalSystem the primal system
             */
            virtual void initializePrimalSystem (dcp::TimeDependentEquationSystem& primalSystem);

            //! Initialize adjoint problem for solving
            /*!
             *  This method relies on \c dcp::TimeDependentProblem::restoreState() method, and it assumes that the
             *  initial state to be restored is called \c "initial_state"
             *
             *  \param adjointSystem the adjoint system
             *  \param primalSystem the primal system, which is not used here but may be useful for derived classes
             */
            virtual void initializeAdjointSystem (dcp::TimeDependentEquationSystem& adjointSystem,
                                                  const dcp::TimeDependentEquationSystem& primalSystem);

            //! Solve the primal problem
            /*!
             *  In this base implementation, \c solve() method is simply called on primalSystem
             *
             *  \param primalSystem the primal system
             */
            virtual void solvePrimalSystem (dcp::TimeDependentEquationSystem& primalSystem);

            //! Solve the adjoint problem
            /*!
             *  It replicates \c dcp::TimeDependentEquationSystem::solve() method but adds the call to
             *  \c linkAdjointToPrimal_() before the system is actually solved at each time step.
             *
             *  \param adjointSystem the adjoint system to be solved
             *  \param primalSystem the primal system, which is then passed to \c linkAdjointToPrimal_()
             */
            virtual void solveAdjointSystem (dcp::TimeDependentEquationSystem& adjointSystem,
                                             const dcp::TimeDependentEquationSystem& primalSystem);

            //! Function to be called before the actual time advancement is performed in \c solvePrimalSystem_()
            /*!
             *  This base method actually does nothing, but it allows derived classes to override just this method if
             *  the problem to be solved is simple enough.
             *
             *  \param primalSystem the primal system
             *  \param adjointSystem the adjoint system
             *  \param timestep the current timestep in the adjoint-system solution loop
             */
            virtual void primalSolveStepPreprocessing (dcp::TimeDependentEquationSystem& primalSystem,
                                                       const dcp::TimeDependentEquationSystem& adjointSystem,
                                                       std::size_t timestep);

            //! Function to be called before the actual time advancement is performed in \c solveAdjointSystem_()
            /*!
             *  This base method actually does nothing, but it allows derived classes to override just this method if
             *  the problem to be solved is simple enough. It should be used, for example, to set the primal
             *  coefficients in the adjoint system
             *
             *  \param adjointSystem the adjoint system
             *  \param primalSystem the primal system
             *  \param timestep the current timestep in the adjoint-system solution loop
             */
            virtual void adjointSolveStepPreprocessing (dcp::TimeDependentEquationSystem& adjointSystem,
                                                        const dcp::TimeDependentEquationSystem& primalSystem,
                                                        std::size_t timestep);

            // ---------------------------------------------------------------------------------------------//

        protected:

            // ---------------------------------------------------------------------------------------------//

        private:

    };



    // ============================================================================================== //
    // ==================================== IMPLEMENTATION ========================================== //
    // ============================================================================================== //

    /******************* CONSTRUCTORS *******************/
    template <class T_ControlVariable>
        TimeDependentBacktrackingImplementer<T_ControlVariable>::TimeDependentBacktrackingImplementer
                (const typename dcp::GenericImplementer<T_ControlVariable>::Updater& updater) :
            GenericImplementer<T_ControlVariable> (updater)
        {
            dolfin::log (dolfin::DBG, "TimeDependentBacktrackingImplementer object created");
        }



    template <class T_ControlVariable>
        TimeDependentBacktrackingImplementer<T_ControlVariable>::TimeDependentBacktrackingImplementer
                (const typename dcp::GenericImplementer<T_ControlVariable>::Updater& updater,
                 const typename dcp::GenericImplementer<T_ControlVariable>::SearchDirectionComputer&
                        searchDirectionComputer) :
            GenericImplementer<T_ControlVariable> (updater, searchDirectionComputer)
        {
            dolfin::log (dolfin::DBG, "TimeDependentBacktrackingImplementer object created");
        }



    /********************** METHODS ***********************/
    template <class T_ControlVariable>
        void TimeDependentBacktrackingImplementer<T_ControlVariable>::solve
            (const std::vector<std::shared_ptr<dcp::GenericEquationSystem> > systems,
             const std::string& solveType)
        {
            // get primal and adjoint problems
            // use of static_cast is justified by the fact that the documentation requires the systems to be of this
            // type
            dcp::TimeDependentEquationSystem& primalSystem =
                static_cast<dcp::TimeDependentEquationSystem&> (*(systems[0]));
            dcp::TimeDependentEquationSystem& adjointSystem =
                static_cast<dcp::TimeDependentEquationSystem&> (*(systems[1]));

            if (solveType != "all" && solveType != "primal" && solveType != "adjoint")
            {
                dolfin::dolfin_error ("dcp: TimeDependentBacktrackingImplementer.cpp",
                                      "solve",
                                      "Unknown solve type \"%s\"",
                                      solveType.c_str ());
            }

            // 1) solve primal problem
            if (solveType == "all" || solveType == "primal")
            {
                dolfin::begin (dolfin::PROGRESS, "Initializing primal system...");
                initializePrimalSystem (primalSystem);
                dolfin::end (); // Initializing primal system

                dolfin::begin (dolfin::PROGRESS, "Solving primal system...");
                solvePrimalSystem (primalSystem);
                dolfin::end (); // Solving primal system
            }

            // 2) solve adjoint problem
            if (solveType == "all" || solveType == "adjoint")
            {
                dolfin::begin (dolfin::PROGRESS, "Initializing adjoint system...");
                initializeAdjointSystem (adjointSystem, primalSystem);
                dolfin::end (); // Initializing adjoint system

                dolfin::begin (dolfin::PROGRESS, "Solving adjoint system...");
                solveAdjointSystem (adjointSystem, primalSystem);
                dolfin::end (); // Solving adjoint system
            }
        }



    template <class T_ControlVariable>
        void TimeDependentBacktrackingImplementer<T_ControlVariable>::initializePrimalSystem
            (dcp::TimeDependentEquationSystem& primalSystem)
        {

            for (std::size_t i = 0; i < primalSystem.size (); ++i)
            {
                dolfin::begin (dolfin::DBG, "Initializing problem number %d in primal problem...", i);
                primalSystem[i].clear ();
                primalSystem[i].restoreState ("initial_state", true);
                dolfin::end (); // Initializing problem number %d in primal problem
            }
        }



    template <class T_ControlVariable>
        void TimeDependentBacktrackingImplementer<T_ControlVariable>::initializeAdjointSystem
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



    template <class T_ControlVariable>
        void TimeDependentBacktrackingImplementer<T_ControlVariable>::solvePrimalSystem
            (dcp::TimeDependentEquationSystem& primalSystem)
        {
            primalSystem.solve ();
        }



    template <class T_ControlVariable>
        void TimeDependentBacktrackingImplementer<T_ControlVariable>::solveAdjointSystem
            (dcp::TimeDependentEquationSystem& adjointSystem,
             const dcp::TimeDependentEquationSystem& primalSystem)
        {
            // i) Reserve space for solutions vector of each problem
            for (std::size_t i = 0; i < adjointSystem.size (); ++i)
            {
                adjointSystem[i].reserve ();
            }

            // ii) Solutions loop
            dolfin::begin (dolfin::DBG, "Solution loop...");
            std::size_t timestep = 0;
            while (adjointSystem.isFinished () == 0)
            {
                timestep++;
                dolfin::begin (dolfin::PROGRESS, "===== Timestep %d =====", timestep);

                adjointSolveStepPreprocessing (adjointSystem, primalSystem, timestep);

                adjointSystem.solve ("step");

                // plot and write to file problems solutions
                for (const auto& problemName : adjointSystem.problemsNames ())
                {
                    // plot solution
                    dcp::TimeDependentProblem& problem = adjointSystem[problemName];
                    int plotInterval = problem.parameters["plot_interval"];

                    if (plotInterval > 0 && timestep % plotInterval == 0)
                    {
                        problem.plotSolution ("last");
                    }

                    // write solution to file
                    int writeInterval = problem.parameters ["write_interval"];
                    if (writeInterval > 0 && timestep % writeInterval == 0)
                    {
                        problem.writeSolutionToFile ("last");
                    }
                }

                dolfin::end (); // Timestep %d
            }

            dolfin::end (); // Solution loop
        }



    template <class T_ControlVariable>
        void TimeDependentBacktrackingImplementer<T_ControlVariable>::primalSolveStepPreprocessing
            (dcp::TimeDependentEquationSystem& primalSystem,
             const dcp::TimeDependentEquationSystem& adjointSystem,
             std::size_t timestep)
        {
        }



    template <class T_ControlVariable>
        void TimeDependentBacktrackingImplementer<T_ControlVariable>::adjointSolveStepPreprocessing
            (dcp::TimeDependentEquationSystem& adjointSystem,
             const dcp::TimeDependentEquationSystem& primalSystem,
             std::size_t timestep)
        {
        }
}

#endif


