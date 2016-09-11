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

#ifndef SRC_OPTIMIZERS_TIMEDEPENDENTBACKTRACKINGOPTIMIZER_H_INCLUDE_GUARD
#define SRC_OPTIMIZERS_TIMEDEPENDENTBACKTRACKINGOPTIMIZER_H_INCLUDE_GUARD

#include <dcp/optimizers/BacktrackingOptimizer.h>
#include <dcp/problems/TimeDependentEquationSystem.h>
#include <dolfin/function/GenericFunction.h>
#include <dolfin/function/Function.h>
#include <functional>
#include <string>
#include <fstream>

namespace dcp
{
    /*! \class TimeDependentBacktrackingOptimizer TimeDependentBacktrackingOptimizer.h
     *  \brief Class that implements the gradient method with backtracking for time dependent optimization. This class
     *  automatically deals with the forward/backward-in-time primal-adjoint system.
     */
    class TimeDependentBacktrackingOptimizer : public dcp::BacktrackingOptimizer
    {
        // ---------------------------------------------------------------------------------------------//

        public:
            /************************* CONSTRUCTORS ********************/
            //! Default constructor
            TimeDependentBacktrackingOptimizer ();

            /************************* DESTRUCTOR ********************/
            //! Destructor
            /*! Default destructor, since members of the class are trivially 
             * destructible.
             */
            virtual ~TimeDependentBacktrackingOptimizer () {};
            
            
            /********************** METHODS ***********************/
            //! Perform optimization on the input problem using the gradient method with backtracking
            /*! 
             *  Input arguments are:
             *  \param systems the systems (possibly more than one) that represent the primal/adjoint system. In this
             *  case, the first two elements of the input vector of systems are supposed to represent the primal 
             *  (forward-in-time) system and the adjoint (backward-in-time) system respectively (the vector may have 
             *  more than one element, but only the first one will be used). The function will set the parameter
             *  \c "purge_interval" to 0 in order to ensure that all the solutions are kept when solving the problems,
             *  since we need all the primal solutions to solve the adjoint system and we may need all the adjoint
             *  solutions to compute the value of the objective functional. 
             *  Both objects must be \c dcp::TimeDependentEquationSystem (or child classes). Initial and final times
             *  must coincide, and the timestep must be equal in magnitude and opposite in sign.
             *  \param objectiveFunctional the objective functional to be minimized
             *  \param initialGuess the starting point for the minimization algorithm. At the end of the function, it
             *  will containt the final value of the control variable
             *  \param updater callable object to update the control parameter value. It can be either be a function 
             *  pointer, a function object or a lambda expression. Its input argument are:
             *  \li the system to update
             *  \li the new value of the control function
             * 
             *  Functors for the most common types of update are provided: see \c dcp::DirichletControlUpdater,
             *  \c dcp::DistributedControlUpdater and \c dcp::NeumannControlUpdater.
             */
            virtual void apply (const std::vector<std::shared_ptr<dcp::GenericEquationSystem> > systems,
                                const dcp::GenericObjectiveFunctional& objectiveFunctional, 
                                dolfin::Function& initialGuess,
                                const dcp::GenericDescentMethod::Updater& updater) override;
            

            // ---------------------------------------------------------------------------------------------//

        protected:
            /********************** METHODS ***********************/
            //! Solve the equation systems representing the primal and the adjoint problem.
            /*!
             *  This function calls the subfunctions \c which have a pretty standard and unrefined behaviour, since
             *  automatizing the treatment of all possible cases would be hard and lead to unintuitive API. If one needs
             *  a specific behaviour, they need to fine tune the class by deriving it and overriding the subfunctions
             *  called by this method.
             *
             *  \param systems the set of systems to be solved
             */
            virtual void solve_ (const std::vector<std::shared_ptr<dcp::GenericEquationSystem> > systems) override;

            //! Initialize primal problem for solving
            /*!
             *  This method relies on \c dcp::TimeDependentProblem::restoreState() method, and it assumes that the
             *  initial state to be restored is called \c "initial_state"
             *
             *  \param primalSystem the primal system
             */
            virtual void initializePrimalSystem_ (dcp::TimeDependentEquationSystem& primalSystem);

            //! Initialize adjoint problem for solving
            /*!
             *  This method relies on \c dcp::TimeDependentProblem::restoreState() method, and it assumes that the
             *  initial state to be restored is called \c "initial_state"
             *
             *  \param adjointSystem the adjoint system
             *  \param primalSystem the primal system, which is not used here but may be useful for derived classes
             */
            virtual void initializeAdjointSystem_ (dcp::TimeDependentEquationSystem& adjointSystem,
                                                   const dcp::TimeDependentEquationSystem& primalSystem);

            //! Solve the primal problem
            /*!
             *  In this base implementation, \c solve() method is simply called on primalSystem
             *
             *  \param primalSystem the primal system
             */
            virtual void solvePrimalSystem_ (dcp::TimeDependentEquationSystem& primalSystem);

            //! Solve the adjoint problem
            /*!
             *  It replicates \c dcp::TimeDependentEquationSystem::solve() method but adds the call to 
             *  \c linkAdjointToPrimal_() before the system is actually solved at each time step.
             *
             *  \param adjointSystem the adjoint system to be solved
             *  \param primalSystem the primal system, which is then passed to \c linkAdjointToPrimal_()
             */
            virtual void solveAdjointSystem_ (dcp::TimeDependentEquationSystem& adjointSystem,
                                              const dcp::TimeDependentEquationSystem& primalSystem);

            //! Function to be called before the actual time advancement is performed in \c solvePrimalSystem_()
            /*!
             *  This base method actually does nothing, but it allows derived classes to override just this method if
             *  the problem to be solved is simple enough.
             *
             *  \param primalSystem the primal system
             *  \param adjointSystem the adjoint system
             *  \param timeStep the current timestep in the adjoint-system solution loop
             */
            virtual void primalSolveStepPreprocessing_ (dcp::TimeDependentEquationSystem& primalSystem,
                                                        const dcp::TimeDependentEquationSystem& adjointSystem,
                                                        unsigned int timeStep);

            //! Function to be called before the actual time advancement is performed in \c solveAdjointSystem_()
            /*!
             *  This base method actually does nothing, but it allows derived classes to override just this method if
             *  the problem to be solved is simple enough. It should be used, for example, to set the primal
             *  coefficients in the adjoint system
             *
             *  \param adjointSystem the adjoint system
             *  \param primalSystem the primal system
             *  \param timeStep the current timestep in the adjoint-system solution loop
             */
            virtual void adjointSolveStepPreprocessing_ (dcp::TimeDependentEquationSystem& adjointSystem,
                                                         const dcp::TimeDependentEquationSystem& primalSystem,
                                                         unsigned int timeStep);

            //! Update the equations systems represeting the primal and the adjoint problem by using the \c updater
            /*!
             *  \param systems the set of systems to be solved
             *  \param updater the functional to be used to update the system (see \c apply() method documentation)
             *  \param control the current value of the control function
             */
            virtual void update_ (const std::vector<std::shared_ptr<dcp::GenericEquationSystem> > systems,
                                  const dcp::GenericDescentMethod::Updater& updater,
                                  const dolfin::GenericFunction& control) override;


            /********************** VARIABLES ***********************/


            // ---------------------------------------------------------------------------------------------//

        private:

    };
}

#endif


