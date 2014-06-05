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

#ifndef SRC_OPTIMIZER_BACKTRACKINGOPTIMIZER_HPP_INCLUDE_GUARD
#define SRC_OPTIMIZER_BACKTRACKINGOPTIMIZER_HPP_INCLUDE_GUARD

#include <Optimizer/AbstractOptimizer.hpp>
#include <ObjectiveFunctional/AbstractObjectiveFunctional.hpp>
#include <dolfin/function/GenericFunction.h>
#include <dolfin/function/Function.h>
#include <functional>
#include <string>

namespace DCP
{
    /*! \class BacktrackingOptimizer BacktrackingOptimizer.hpp
     *  \brief Class that implements the gradient method with backtracking.
     *  
     *  This class provides the implementation of the gradient method as a descent method
     *  for optimazation of funcionals. It is derived from \c DCP::AbstractOptimizer.
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
    class BacktrackingOptimizer : public DCP::AbstractOptimizer
    {
        // ---------------------------------------------------------------------------------------------//
        
        public:
            /************************* CONSTRUCTORS ********************/
            //! Default constructor
            /*!
             *  Input argument:
             *  \param gradientNormTolerance the tolerance for convergence check on gradient (see third input parameter).
             *  Default value: \c 1e-6
             *  \param relativeIncrementTolerance the tolerance for convergence check on increment (see third input parameter)
             *  Default value: \c 1e-6
             *  \param convergenceCriterion 
             *  sets the convergence criterion to stop the minimization algorithm.
             *  Possible values are:
             *  \li \c gradient: the minimization loop permanence condition is
             *  \f$ \left| \left| \nabla J \right| \right| > \varepsilon_g \f$
             *  where \f$ \varepsilon_g \f$ is set by the input variable \c gradientNormTolerance
             *  \li \c increment: the minimization loop permanence condition is
             *  \f$ \frac{\left| \left| \mathbf{u}_{k+1} - \mathbf{u}_k \right| \right|}
             *           {\left| \left| \mathbf{u}_{k} \right| \right|} 
             *           > \varepsilon_i \f$
             *  where \f$ \varepsilon_i \f$ is set by the input variable \c relativeIncrementTolerance
             *  \li \c both: both the above conditions are checked for the permanence in the minimization loop. That means
             *  that the minimization algorithm will end when one of the two conditions is false
             *  
             *  Default value: \c both
             */
            BacktrackingOptimizer (const double& gradientNormTolerance = 1e-6,
                                   const double& relativeIncrementTolerance = 1e-6,
                                   const std::string& convergenceCriterion = "both");
            
            /************************* DESTRUCTOR ********************/
            //! Destructor
            /*! Default destructor, since members of the class are trivially 
             * destructible.
             */
            virtual ~BacktrackingOptimizer () {};
            
            
            /********************** METHODS ***********************/
            //! Sets the form used to compute dot products and norms. 
            /*!
             *  \param dotProductComputer the form to be stored in the protected member \c dotProductComputer_
             */
            virtual void setDotProductComputer (const boost::shared_ptr<dolfin::Form> dotProductComputer);
            
            //! Reset the value of the protected membet \c dotProductComputer_ so that the default form will be used
            virtual void resetDotProductComputer ();
                
            //! Perform optimization on the input problem using the gradient method with backtracking
            /*! 
             *  Input arguments are:
             *  \param problem the composite differential problem that represents the primal/adjoint system
             *  \param objectiveFunctional the objective functional to be minimized
             *  \param initialGuess the starting point for the minimization algorithm. At the end of the function, it
             *  will containt the final value of the control variable
             *  \param updater callable object to update the control parameter value. It can be either be a function 
             *  pointer, a function object or a lambda expression. Its input argument are:
             *  \li the composite differential problem to update
             *  \li the new value of the control function
             *  
             *  \param searchDirectionComputer callable object to compute the search direction. It can either be a 
             *  function pointer, a function object or a lambda expression. In general the search direction is computed
             *  as: 
             *  \f[
             *      \mathbf{d}_k = -B_k\,\nabla J_k
             *  \f]
             *  The default value is the member function \c gradientSearchDirection(), that basically uses the above 
             *  formula with \f$ B_k = I \f$.
             *  The input arguments for \c searchDirectionComputer are:
             *  \li the dolfin function that will contain the search direction after the function exits
             *  \li the dolfin function containing the gradient
             */
            virtual void apply (DCP::CompositeDifferentialProblem& problem,
                                const DCP::AbstractObjectiveFunctional& objectiveFunctional, 
                                dolfin::Function& initialGuess,
                                const std::function 
                                <
                                    void (DCP::CompositeDifferentialProblem&, const dolfin::GenericFunction&)
                                >& updater,
                                const std::function
                                <
                                    void (dolfin::Function&, const dolfin::Function&)
                                >& searchDirectionComputer = BacktrackingOptimizer::gradientSearchDirection);
            
            //! Perform optimization on the input problem using the gradient method with backtracking and dumping
            //! results to file in the process
            /*! 
             *  Input arguments are:
             *  \param problem the composite differential problem that represents the primal/adjoint system
             *  \param objectiveFunctional the objective functional to be minimized
             *  \param initialGuess the starting point for the minimization algorithm. At the end of the function, it
             *  will containt the final value of the control variable
             *  \param updater callable object to update the control parameter value. It can be either be a function 
             *  pointer, a function object or a lambda expression. Its input argument are:
             *  \li the composite differential problem to update
             *  \li the new value of the control function
             *  
             *  \param dumper callable object to peform the dump of data during the minimization iterations. It can 
             *  either be a function pointer, a function object or a lambda expression. 
             *  \param dumpInterval integer that defines the frequency of dumping. The \c dumper will be called
             *  every \c dumpInterval minimization iterations (and after the last iteration). 
             *  If \c dumpInterval is set to \c 0, the \c dumper will be called only after the last iteration.
             *  \param searchDirectionComputer callable object to compute the search direction. It can either be a 
             *  function pointer, a function object or a lambda expression. In general the search direction is computed
             *  as: 
             *  \f[
             *      \mathbf{d}_k = -B_k\,\nabla J_k
             *  \f]
             *  The default value is the member function \c gradientSearchDirection(), that basically uses the above 
             *  formula with \f$ B_k = I \f$.
             *  The input arguments for \c searchDirectionComputer are:
             *  \li the dolfin function that will contain the search direction after the function exits
             *  \li the dolfin function containing the gradient
             */
            virtual void apply (DCP::CompositeDifferentialProblem& problem,
                                const DCP::AbstractObjectiveFunctional& objectiveFunctional, 
                                dolfin::Function& initialGuess,
                                const std::function 
                                <
                                    void (DCP::CompositeDifferentialProblem&, const dolfin::GenericFunction&)
                                >& updater,
                                const std::function 
                                <
                                    void ()
                                >& dumper,
                                const int& dumpInterval,
                                const std::function
                                <
                                    void (dolfin::Function&, const dolfin::Function&)
                                >& searchDirectionComputer = BacktrackingOptimizer::gradientSearchDirection);
            
            //! Function to compute the search direction. See documentation of method \c apply() for more information
            static void gradientSearchDirection (dolfin::Function& searchDirection, const dolfin::Function& gradient);

            // ---------------------------------------------------------------------------------------------//

        protected:
            //! Function to print to file some values. It is mostly useful to avid code repetition inside this class
            void print (std::ostream& OUTSTREAM, 
                        const int& iteration,
                        const double& functionalValue,
                        const double& alpha,
                        const int& backtrackingIterations,
                        const double& gradientNorm,
                        const double& relativeIncrement);
            
            //! Empty dumper, used for compatibility. This function will be used when the version that does not 
            //! take a \c dumper as input argument of the method \c apply() is called
            static void emptyDumper () 
            {
                dolfin::log (dolfin::DBG, "Empty dumper was called");
            };
                
            //! The form that will be used to compute the dot product between the gradient and the search direction. 
            /*! 
             *  The default value is \c nullptr, in which case the function will use one
             *  of the forms in \c DotProduct.h, trying to determine the right one by checking the geometrical dimensions
             *  of the input objects. However, sometimes it may be useful to pass a user-defined object to perform
             *  the task. To do so, use the function \c setDotProductComputer
             */
            boost::shared_ptr<dolfin::Form> dotProductComputer_;

            // ---------------------------------------------------------------------------------------------//

        private:
            
    };
}

#endif

