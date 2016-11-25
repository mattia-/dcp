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

#ifndef SRC_OPTIMIZERS_BACKTRACKINGOPTIMIZER_H_INCLUDE_GUARD
#define SRC_OPTIMIZERS_BACKTRACKINGOPTIMIZER_H_INCLUDE_GUARD

#include <functional>
#include <string>
#include <fstream>
#include <dolfin/function/GenericFunction.h>
#include <dolfin/function/Function.h>
#include <dcp/optimizers/GenericDescentMethod.h>
#include <dcp/objective_functionals/GenericObjectiveFunctional.h>
#include <dcp/optimizers/GenericImplementer.h>

namespace dcp
{
    /*! \class BacktrackingOptimizer BacktrackingOptimizer.h
     *  \brief Class that implements the gradient method with backtracking.
     *
     *  This class provides the implementation of the gradient method as a descent method
     *  for optimazation of functionals. It is derived from \c dcp::GenericDescentMethod.
     *  Let \f$ J \f$ be the functional to be minimized and \f$ \psi \left( \alpha \right) \f$ the function defined as:
     *  \f[
     *      \psi \left( \alpha \right) = J \left( \mathbf{u} + \alpha\,\mathbf{d} \right)
     *  \f]
     *  with \f$ \mathbf{d} \f$ search direction.
     *  As for all descent method, the generic iteration of the gradient method with backtracking can be
     *  written as:
     *  \f[
     *      \mathbf{u}_{k+1} = \mathbf{u}_k + \alpha_k\,\mathbf{d}_k\,.
     *  \f]
     *  The search direction \f$\mathbf{d}_k\f$ is chosen to be:
     *  \f[
     *      \mathbf{d}_k = - \nabla J \left(\mathbf{u}_k\right)
     *  \f]
     *  The gradient method with backtracking performs an approximated line minimization to determine the
     *  weigth \f$ \alpha \f$, imposing only the first Wolfe condition (i.e. the sufficient-decrease condition,
     *  aka Armijo Rule) and using a backtracking algorithm to determine \f$ \alpha_k \f$, the weight of the
     *  gradient.
     *  In particular, the algorithm is as follows. \n
     *  UNTIL convergence, for every \f$ k \f$: \n
     *  let \f$ c_1 \f$, \f$ \alpha^{\left(0\right)}_k \f$ and \f$ \rho \f$ be given and set \f$ m = 0 \f$. Then: \n
     *  WHILE \f$ \left( \psi \left( \alpha ^ {\left( m \right)}_k \right) >
     *                   \psi \left(0\right) + c_1\,\alpha^{\left( m \right)}_k\,\psi' \left(0\right)
     *            \right) \f$
     *  \li \f$ \alpha^{\left(m + 1\right)}_k = \rho \alpha^{\left(m\right)}_k \f$ \n
     *  \li \f$ m++ \f$
     *
     *  All the parameters are stored in the class member \c parameters and are given default values:
     *  \li \f$ c_1 = 10^{-3}\f$
     *  \li \f$ \alpha^{\left(0\right)} = 0.5 \f$
     *  \li \f$ \rho = 0.5 \f$
     */
    class BacktrackingOptimizer : public dcp::GenericDescentMethod
    {
        // ---------------------------------------------------------------------------------------------//

        public:
            /************************* CONSTRUCTORS ********************/
            //! Default constructor
            /*!
             *  The constructors also sets the following parameters:
             *      - \c "descent_method" a name identifying this algorithm. Default value: \c backtracking_gradient_method
             *      - \c "gradient_norm_tolerance" the tolerance for convergence check on gradient. Default value: \c 1e-6
             *      - \c "relative_increment_tolerance" the tolerance for convergence check on increment.
             *        Default value: \c 1e-6
             *      - \c "convergence_criterion" sets the convergence criterion to stop the minimization algorithm.
             *        Possible values are:
             *           - \c \c "gradient": the minimization loop permanence condition is
             *             \f$ \left| \left| \nabla J \right| \right| > \varepsilon_g \f$
             *             where \f$ \varepsilon_g \f$ is set by the parameter \c gradient_norm_tolerance
             *           - \c \c "increment": the minimization loop permanence condition is
             *             \f$ \frac{\left| \left| \mathbf{u}_{k+1} - \mathbf{u}_k \right| \right|}
             *                       {\left| \left| \mathbf{u}_{k} \right| \right|}
             *                       > \varepsilon_i \f$
             *             where \f$ \varepsilon_i \f$ is set by the parameter \c relative_increment_tolerance
             *           - \c \c "both": both the above conditions are checked for the permanence in the minimization
             *             loop. That means that the minimization algorithm will end when one of the two conditions is
             *             false.
             *        Default value: \c both
             *      - \c "c_1" the value of the parameter \c c_1 in the backtracking algorithm (see class documentation).
             *        Default value: 1e-3
             *      - \c "alpha" the value of the parameter \c alpha in the backtracking algorithm (see class
             *        documentation). Default value: 1.0
             *      - \c "rho" the value of the parameter \c rho in the backtracking algorithm (see class documentation).
             *        Default value: 0.5
             *      - \c "max_minimization_iterations" the maximum number of iteration for the minimization loop.
             *        Default value: 100
             *      - \c "max_backtracking_iterations" the maximum number of iteration for the backtracking loop.
             *        Default value: 20
             *      - \c "output_file_name" the name of the output file to write results to.
             *        Default value: empty (which means that output will be written to terminal only).
             */
            BacktrackingOptimizer ();

            /************************* DESTRUCTOR ********************/
            //! Destructor
            /*! Default destructor, since members of the class are trivially
             * destructible.
             */
            virtual ~BacktrackingOptimizer () {};


            /********************** METHODS ***********************/
            //! Perform optimization on the input problem using the gradient method with backtracking
            /*!
             *  Just a wrapper to allow easier calls to the real \c apply() method, which is the one that takes a vector
             *  of systems as input. See that function's documentation for in-depth explaination of the arguments
             *
             *  Input arguments are:
             *  \param system the system that represent the primal/adjoint system
             *  \param objectiveFunctional the objective functional to be minimized
             *  \param initialGuess the starting point for the minimization algorithm. At the end of the function, it
             *  will contain the final value of the control variable
             *  \param implementer object that defines all the methods needed by the algorithm
             */
            template <class T_ControlVariable>
                void apply (dcp::GenericEquationSystem& system,
                            dcp::GenericObjectiveFunctional& objectiveFunctional,
                            T_ControlVariable& initialGuess,
                            dcp::GenericImplementer<T_ControlVariable>& implementer) const;

            //! Perform optimization on the input problem using the gradient method with backtracking
            /*!
             *  The type \c T_ControlVariable must satisfy the following constraints:
             *  \li be copy-constructible
             *  \li have a sum operator and a multiplication by scalar operator defined
             *
             *
             *  Input arguments are:
             *  \param systems the systems (possibly more than one) that represent the primal/adjoint system. There may
             *  be more than one for example in time-dependent cases, where the primal and adjoint systems must be
             *  implemented separately since they have opposite direction of time advancement; if the standard
             *  implementers are used (see documentation for \c implementer ) the primal-adjoint system is supposed to
             *  be stored in the first element of the vector in the stationary case and in the first two elements
             *  (primal in the first and adjoint in the second) of the vector in the time dependent case
             *  \param objectiveFunctional the objective functional to be minimized
             *  \param initialGuess the starting point for the minimization algorithm. At the end of the function, it
             *  will contain the final value of the control variable
             *  \param implementer object that defines all the methods needed by the algorithm. The
             *  \c apply() method is indeed general and delegates the implementations of the particular methods to the
             *  \c implementer . Some implementers are already given as part of \c dcp ,  namely 
             *  \c dcp::BacktrackingImplementer and \c dcp::TimeDependentBacktrackingImplementer . The object passed to
             *  the \c apply() function may be an instance of either these classes or a class derived from them or even
             *  a completely new class. The class of which this object is an instance must implement the following 
             *  methods:
             *  \li <tt>void update (const std::vector<std::shared_ptr<dcp::GenericEquationSystem>, 
             *                       const T_ControlVariable&)</tt>
             *  method that updates the systems (passed as first argument) using the new value of the control function
             *  passed as second argument; it is called after each time the control variable changes value. Some
             *  standard functors that can be used within the implementer are provided (see 
             *  \c dcp::DirichletControlUpdater , \c dcp::DistributedControlUpdater and \c dcp::NeumannControlUpdater )
             *  \li <tt>T_ControlVariable computeSearchDirection (const T_ControlVariable&)</tt>
             *  method that computes the search direction given the functional gradient (passed as first and only
             *  argument); note that C++11 automatically implements the move semantic for object returned from functions
             *  so there is no overhead even in returning a big object (hopefully); it is called at the beginning of
             *  each minimization iteration
             *  \li <tt>void solve (const std::vector<std::shared_ptr<dcp::GenericEquationSystem> >,
             *                      const std::string&)</tt>
             *  method that solves the problems stored in the argument \c systems (passed as first argument to the
             *  function); the second argument contains a string switch to (possiby) choose between different behaviours
             *  \li <tt>double computeDotProduct (const T_ControlVariable&, const T_ControlVariable&)</tt>
             *  method to compute the dot product between to functions of type \c T_ControlVariable ; it is called on
             *  the functional gradient (first argument) and the search direction (second argument), and its result will
             *  be used in the sufficient decrease condition.
             *  \li <tt>double computeNorm (const T_ControlVariable&)</tt>
             *  method to compute the norm of the gradient or of the increment of the control variable (passed as 
             *  argument to the function); it is used for the convergence check, if needed. Note that this could use 
             *  the previous \c computeDotProduct function within the implementer
             */
            template <class T_ControlVariable>
                void apply (const std::vector<std::shared_ptr<dcp::GenericEquationSystem> > systems,
                            dcp::GenericObjectiveFunctional& objectiveFunctional,
                            T_ControlVariable& initialGuess,
                            dcp::GenericImplementer<T_ControlVariable>& implementer) const;

            // ---------------------------------------------------------------------------------------------//

        protected:
            //! Function to perform the backtracking loop iterations
            /*
             *  \param previousFunctionalValue the value of the functional at the beginning of the backtracking loop
             *  \param currentFunctionalValue the new value of the functional on each iteration. At the end of the
             *  function, it will contain the value of the functional at the last iteration.
             *  \param gradientDotSearchDirection the dot product between the functional gradient and the search
             *  direction. Used for the sufficient decrease condition.
             *  \param alpha the value of alpha in the backtracking loop. At the end of the function, it will contain
             *  the value of the alpha at the last iteration.
             *  \param backtrackingIteration the iterations performed in the backtracking loop. At the end of the
             *  function, it will contain the total number of iterations performed.
             *  \param controlVariable the control variable. At the end of the function, it will contain the
             *  control variable at the last iteration.
             *  \param previousFunctionalValue the control variable at the beginning of the backtracking loop.
             *  \param searchDirection the search direction
             *  \param systems the systems that represent the primal/adjoint system
             *  \param objectiveFunctional the objective functional
             *  \param implementer the implementer (see the documentation for the method \c apply() for a deeper
             *  explanation)
             *
             *  \return \c true if the sufficient decrease condition is satisfied at the exit of the loop,
             *  \c false otherwise
             */
            template <class T_ControlVariable>
            bool backtrackingLoop_ (const double& previousFunctionalValue,
                                    double& currentFunctionalValue,
                                    const double& gradientDotSearchDirection,
                                    double& alpha,
                                    int& backtrackingIteration,
                                    T_ControlVariable& controlVariable,
                                    const T_ControlVariable& previousControlVariable,
                                    const T_ControlVariable& searchDirection,
                                    const std::vector<std::shared_ptr<dcp::GenericEquationSystem> > systems,
                                    dcp::GenericObjectiveFunctional& objectiveFunctional,
                                    dcp::GenericImplementer<T_ControlVariable>& implementer) const;

            //! Function to print to file some values. It is mostly useful to avoid code repetition inside this class
            /*!
             *  Input parameters:
             *  \param iteration the current iteration when the function is called
             *  \param OUTSTREAM the stream to print to
             *  \param functionalValue the value of the functional
             *  \param alpha the backtracking parameter
             *  \param backtrackingIterations the number of iterations performed in the backtracking loop
             *  \param gradientNorm the norm of the gradient of the functional
             *  \param relativeIncrement the relative increment of the optimal solution found
             */
            void print_ (std::ostream& OUTSTREAM,
                         const int& iteration,
                         const double& functionalValue,
                         const double& alpha,
                         const int& backtrackingIterations,
                         const double& gradientNorm,
                         const double& relativeIncrement) const;

            //! Function to open output file and print the header
            /*!
             *  \param outfile the file stream to print to
             *
             *  \return \c true if the file was opened, \c false otherwise
             */
            bool openOutputFile_ (std::ofstream& outfile) const;

            /********************** VARIABLES ***********************/


            // ---------------------------------------------------------------------------------------------//

        private:

    };



    // ============================================================================================== //
    // ==================================== IMPLEMENTATION ========================================== //
    // ============================================================================================== //

    template <class T_ControlVariable>
        void BacktrackingOptimizer::apply (dcp::GenericEquationSystem& system,
                                           dcp::GenericObjectiveFunctional& objectiveFunctional,
                                           T_ControlVariable& initialGuess,
                                           dcp::GenericImplementer<T_ControlVariable>& implementer) const
        {
            std::vector<std::shared_ptr<dcp::GenericEquationSystem> > systemAsVector;
            systemAsVector.push_back (dolfin::reference_to_no_delete_pointer (system));

            apply (systemAsVector, objectiveFunctional, initialGuess, implementer);
        }



    template <class T_ControlVariable>
        void BacktrackingOptimizer::apply (const std::vector<std::shared_ptr<dcp::GenericEquationSystem> > systems,
                                           dcp::GenericObjectiveFunctional& objectiveFunctional,
                                           T_ControlVariable& initialGuess,
                                           dcp::GenericImplementer<T_ControlVariable>& implementer) const
    {
        // ------------------------------------------------------------------------------------------------------- //
        // minimization loop settings
        // ------------------------------------------------------------------------------------------------------- //
        // get parameters values
        double c_1                        = this->parameters ["c_1"];
        double alpha_0                    = this->parameters ["alpha_0"];
        double rho                        = this->parameters ["rho"];
        double gradientNormTolerance      = this->parameters ["gradient_norm_tolerance"];
        double relativeIncrementTolerance = this->parameters ["relative_increment_tolerance"];
        std::string convergenceCriterion  = this->parameters ["convergence_criterion"];
        int maxMinimizationIterations     = this->parameters ["max_minimization_iterations"];
        int maxBacktrackingIterations     = this->parameters ["max_backtracking_iterations"];

        // define loop variables
        double alpha;
        double gradientNorm;
        double gradientDotSearchDirection;
        double incrementNorm;
        double previousControlVariableNorm;
        double relativeIncrement;
        double previousFunctionalValue;
        double currentFunctionalValue;
        T_ControlVariable& controlVariable = initialGuess;

        // the following initializations have no purpose but to create variables of the same type as initialGuess
        T_ControlVariable functionalGradient = initialGuess;
        T_ControlVariable searchDirection = initialGuess;
        T_ControlVariable controlVariableIncrement = initialGuess;

        // define output file and print header if necessary
        std::ofstream outfile;
        bool hasOutputFile = openOutputFile_ (outfile);

        // all linear problems in system[0] should be reassembled every time. So we set the parameter
        // "force_reassemble_system" to true for every one of them
        dolfin::begin (dolfin::DBG, "Setting parameter \"force_reassemble_system\" to true for all problems...");
        for (std::size_t i = 0; i < systems.size (); ++i)
        {
            for (std::size_t j = 0; j < systems[i]->size (); ++j)
            {
                if ((*(systems[i]))[j].parameters.has_key ("force_reassemble_system"))
                {
                    (*(systems[i]))[j].parameters["force_reassemble_system"] = true;
                }
            }
        }
        dolfin::end ();

        // set purge_inteval equal to 0, because we need all the solutions when we solve the backward-in-time adjoint
        // system 
        dolfin::begin (dolfin::DBG, "Setting parameter \"purge_interval\" to 0 for all problems...");
        for (std::size_t i = 0; i < systems.size (); ++i)
        {
            for (std::size_t j = 0; j < systems[i]->size (); ++j)
            {
                if ((*(systems[i]))[j].parameters.has_key ("purge_interval"))
                {
                    (*(systems[i]))[j].parameters["purge_interval"] = 0;
                }
            }
        }
        dolfin::end ();

        // function to check convergence
        std::function <bool ()> isConverged;

        // set value of isConverged () according to convergence criterion
        if (convergenceCriterion == "gradient")
        {
            isConverged = [&] () {return gradientNorm < gradientNormTolerance;};
        }
        else if (convergenceCriterion == "increment")
        {
            isConverged = [&] () {return relativeIncrement < relativeIncrementTolerance;};
        }
        else if (convergenceCriterion == "both")
        {
            isConverged = [&] ()
            {
                return (relativeIncrement < relativeIncrementTolerance) || (gradientNorm < gradientNormTolerance);
            };
        }
        else
        {
            dolfin::dolfin_error ("dcp: BacktrackingOptimizer.cpp",
                                  "apply",
                                  "Unknown convergence criterion \"%s\"",
                                  convergenceCriterion.c_str ());
        }
        // ------------------------------------------------------------------------------------------------------- //
        // end of minimization loop setting
        // ------------------------------------------------------------------------------------------------------- //


        // update and solve the system for the first time
        dolfin::begin (dolfin::PROGRESS, "Minimization loop initialization...");

        dolfin::begin (dolfin::PROGRESS, "Updating systems using control initial guess...");
        implementer.update (systems, controlVariable);
        dolfin::end (); // Updating systems using control initial guess

        // solve problems
        dolfin::begin (dolfin::PROGRESS, "Solving...");
        implementer.solve (systems, "all");
        dolfin::end (); // Solving

        // initialize loop variables
        dolfin::begin (dolfin::PROGRESS, "Computing norm of functional gradient...");
        functionalGradient = objectiveFunctional.gradient ();
        gradientNorm = implementer.computeNorm (functionalGradient);
        dolfin::end (); // Computing norm of functional gradient

        relativeIncrement = relativeIncrementTolerance + 1; // just an initialization to be sure that the first iteration
                                                            // of the minimization loop is performed
        dolfin::begin (dolfin::PROGRESS, "Evaluating functional...");
        currentFunctionalValue = objectiveFunctional.evaluateFunctional ();
        dolfin::end (); // Evaluating functional

        dolfin::end (); // "Minimization loop initialization"

        dolfin::log (dolfin::DBG, "==========================\nMINIMIZATION LOOP SETTINGS\n==========================");
        dolfin::log (dolfin::DBG, "Parameters value are:");
        dolfin::log (dolfin::DBG, "c_1 = %f", c_1);
        dolfin::log (dolfin::DBG, "alpha_0 = %f", alpha_0);
        dolfin::log (dolfin::DBG, "rho = %f", rho);
        dolfin::log (dolfin::DBG, "gradient norm tolerance = %f", gradientNormTolerance);
        dolfin::log (dolfin::DBG, "relative increment tolerance = %f", relativeIncrementTolerance);
        dolfin::log (dolfin::DBG, "convergence criterion = %s", convergenceCriterion.c_str ());
        dolfin::log (dolfin::DBG, "maximum minimization iterations = %d", maxMinimizationIterations);
        dolfin::log (dolfin::DBG, "maximum backtracking iterations = %d", maxBacktrackingIterations);


        dolfin::log (dolfin::PROGRESS,
                     "\n============================\nLOOP VARIABLES INITIAL VALUE\n============================");
        dolfin::log (dolfin::PROGRESS, "gradient norm = %f", gradientNorm);
        dolfin::log (dolfin::PROGRESS, "functional value = %f\n", currentFunctionalValue);

        dolfin::log (dolfin::PROGRESS, "***************************");
        dolfin::log (dolfin::PROGRESS, "**** MINIMIZATION LOOP ****");
        dolfin::begin (dolfin::PROGRESS, "***************************");
        int minimizationIteration = 0;

        // print results to file
        if (hasOutputFile)
        {
            print_ (outfile, minimizationIteration, currentFunctionalValue, 0, 0, gradientNorm, relativeIncrement);
        }

        while (isConverged () == false && minimizationIteration < maxMinimizationIterations)
        {
            minimizationIteration++;

            dolfin::log (dolfin::PROGRESS, "==========================");
            dolfin::log (dolfin::PROGRESS, "Minimization iteration %d", minimizationIteration);
            dolfin::begin (dolfin::PROGRESS, "==========================");

            // iteration-specific variable initialization
            alpha = alpha_0;
            int backtrackingIteration = 0;
            previousFunctionalValue = currentFunctionalValue;

            // get current functional gradient;
            // note both primal and adjoint problem have been solver before the loop on the first iteration and then
            // at the end of the previous iteration, so the functional gradient is coherent
            functionalGradient = objectiveFunctional.gradient ();
            dolfin::begin (dolfin::PROGRESS, "Computing search direction...");
            searchDirection = implementer.computeSearchDirection (functionalGradient);
            dolfin::end ();

            // compute dot product between gradient and search direction
            dolfin::begin (dolfin::DBG, "Computing dot product between gradient and search direction...");
            gradientDotSearchDirection = implementer.computeDotProduct (functionalGradient, searchDirection);
            dolfin::end ();

            // ------------------------ solution of system with alpha_0 ------------------------ //
            dolfin::log (dolfin::PROGRESS, "Alpha = %f", alpha);

            // update control variable value
            dolfin::begin (dolfin::PROGRESS, "Updating control variable...");
            T_ControlVariable previousControlVariable (controlVariable);
            controlVariable = previousControlVariable + (searchDirection * alpha);
            dolfin::end ();

            // update system
            dolfin::begin (dolfin::PROGRESS, "Updating systems using control initial guess...");
            implementer.update (systems, controlVariable);
            dolfin::end (); // Updating systems using control initial guess

            // solve system (only primal problem)
            dolfin::begin (dolfin::PROGRESS, "Solving...");
            implementer.solve (systems, "primal");
            dolfin::end (); // Solving

            // update value of the functional
            dolfin::begin (dolfin::PROGRESS, "Evaluating functional...");
            currentFunctionalValue = objectiveFunctional.evaluateFunctional ();
            dolfin::end ();

            dolfin::log (dolfin::PROGRESS, "Functional value = %f\n", currentFunctionalValue);
            // ------------------------ end of solution of system with alpha_0 ------------------------ //

            // backtracking loop
            backtrackingLoop_ (previousFunctionalValue,
                               currentFunctionalValue,
                               gradientDotSearchDirection,
                               alpha,
                               backtrackingIteration,
                               controlVariable,
                               previousControlVariable,
                               searchDirection,
                               systems,
                               objectiveFunctional,
                               implementer);

            // solve adjoint problem after backtracking, so we get the updated functional gradient
            dolfin::begin (dolfin::PROGRESS, "Solving...");
            implementer.solve (systems, "adjoint");
            dolfin::end (); // Solving

            // update gradient norm if necessary for convergence check
            if (convergenceCriterion == "gradient" || convergenceCriterion == "both")
            {
                dolfin::begin (dolfin::DBG, "Computing norm of functional gradient...");
                functionalGradient = objectiveFunctional.gradient ();
                gradientNorm = implementer.computeNorm (functionalGradient);
                dolfin::end ();

                dolfin::log (dolfin::PROGRESS, "Gradient norm = %f", gradientNorm);
            }


            // update relative increment if necessary for convergence check
            if (convergenceCriterion == "increment" || convergenceCriterion == "both")
            {
                dolfin::begin (dolfin::DBG, "Computing relative increment...");
                previousControlVariableNorm = implementer.computeNorm (previousControlVariable);
                dolfin::end ();

                controlVariableIncrement = controlVariable - previousControlVariable;
                incrementNorm = implementer.computeNorm (controlVariableIncrement);

                // compute relative increment. We add DOLFIN_EPS at the denominator in case previousControlVariableNorm
                // is zero (like at the first iteration)
                relativeIncrement = incrementNorm / (previousControlVariableNorm + DOLFIN_EPS);
                dolfin::log (dolfin::PROGRESS, "Relative increment = %f", relativeIncrement);
            }

            dolfin::log (dolfin::PROGRESS, "Functional value = %f\n", currentFunctionalValue);

            if (hasOutputFile)
            {
                dolfin::log (dolfin::DBG, "Printing results to file...");
                print_ (outfile,
                        minimizationIteration,
                        currentFunctionalValue,
                        alpha,
                        backtrackingIteration,
                        gradientNorm,
                        relativeIncrement);
            }

            dolfin::end (); // "Minimization iteration %d"
        }

        dolfin::end (); // "Minimization loop"


        if (isConverged () == false)
        {
            dolfin::log (dolfin::PROGRESS, "End of Minimization loop");
            dolfin::warning ("Maximum number of iterations reached");
            dolfin::log (dolfin::PROGRESS, "Iterations performed: %d", minimizationIteration);
            dolfin::log (dolfin::PROGRESS, "Gradient norm = %f", gradientNorm);
            dolfin::log (dolfin::PROGRESS, "Relative increment = %f", relativeIncrement);
            dolfin::log (dolfin::PROGRESS, "Functional value = %f\n", currentFunctionalValue);
        }
        else
        {
            dolfin::log (dolfin::PROGRESS, "End of Minimization loop");
            dolfin::log (dolfin::PROGRESS, "Iterations performed: %d", minimizationIteration);
            dolfin::log (dolfin::PROGRESS, "Gradient norm = %f", gradientNorm);
            dolfin::log (dolfin::PROGRESS, "Relative increment = %f", relativeIncrement);
            dolfin::log (dolfin::PROGRESS, "Functional value = %f\n", currentFunctionalValue);
        }

        dolfin::log (dolfin::PROGRESS, "-----------------------------\n");

        if (hasOutputFile)
        {
            outfile.close ();
        }
    }



    template <class T_ControlVariable>
        bool BacktrackingOptimizer::backtrackingLoop_
            (const double& previousFunctionalValue,
             double& currentFunctionalValue,
             const double& gradientDotSearchDirection,
             double& alpha,
             int& backtrackingIteration,
             T_ControlVariable& controlVariable,
             const T_ControlVariable& previousControlVariable,
             const T_ControlVariable& searchDirection,
             const std::vector<std::shared_ptr<dcp::GenericEquationSystem> > systems,
             dcp::GenericObjectiveFunctional& objectiveFunctional,
             dcp::GenericImplementer<T_ControlVariable>& implementer) const
    {
        // get parameters
        double c_1 = this->parameters ["c_1"];
        double rho = this->parameters ["rho"];
        int maxBacktrackingIterations = this->parameters ["max_backtracking_iterations"];

        // function to check if sufficient-decrease condition is satisfied
        auto sufficientDecreaseConditionIsSatisfied = [&c_1, &gradientDotSearchDirection]
            (const double& previousValue, const double& currentValue, const double& alpha)
        {
            return currentValue < previousValue + c_1 * alpha * gradientDotSearchDirection;
        };

        dolfin::begin (dolfin::PROGRESS, "Starting backtracking loop...");
        while (sufficientDecreaseConditionIsSatisfied (previousFunctionalValue, currentFunctionalValue, alpha) == false
               && backtrackingIteration < maxBacktrackingIterations)
        {
            // iteration-specific variable initialization
            ++backtrackingIteration;
            alpha = alpha * rho;

            dolfin::log (dolfin::PROGRESS, "========== Backtracking iteration %d ==========", backtrackingIteration);
            dolfin::log (dolfin::PROGRESS, "Alpha = %f", alpha);

            // update control variable value
            dolfin::begin (dolfin::PROGRESS, "Updating control variable...");
            controlVariable = previousControlVariable + (searchDirection * alpha);
            dolfin::end ();

            // update system
            dolfin::begin (dolfin::PROGRESS, "Updating differential problem...");
            implementer.update (systems, controlVariable);
            dolfin::end ();

            // solve system (only primal problem)
            dolfin::begin (dolfin::PROGRESS, "Solving...");
            implementer.solve (systems, "primal");
            dolfin::end (); // Solving

            // update value of the functional
            dolfin::begin (dolfin::PROGRESS, "Evaluating functional...");
            currentFunctionalValue = objectiveFunctional.evaluateFunctional ();
            dolfin::end ();

            dolfin::log (dolfin::PROGRESS, "Functional value = %f", currentFunctionalValue);

            dolfin::log (dolfin::PROGRESS, "");
        }

        dolfin::end ();

        if (sufficientDecreaseConditionIsSatisfied (previousFunctionalValue, currentFunctionalValue, alpha) == false)
        {
            dolfin::log (dolfin::PROGRESS, "");
            dolfin::warning ("Backtracking loop ended because maximum number of iterations was reached");
            dolfin::log (dolfin::PROGRESS, "Alpha (determined with backtracking) = %f", alpha);
            return false;
        }
        else
        {
            dolfin::log (dolfin::PROGRESS, "");
            dolfin::log (dolfin::PROGRESS, "Backtracking loop ended. Iterations performed: %d\n", backtrackingIteration);
            dolfin::log (dolfin::PROGRESS, "Alpha (determined with backtracking) = %f\n", alpha);
            return true;
        }
    }

}

#endif
