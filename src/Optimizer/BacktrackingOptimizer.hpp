#ifndef SRC_DESCENTMETHOD_GRADIENTMETHOD_HPP_INCLUDE_GUARD
#define SRC_DESCENTMETHOD_GRADIENTMETHOD_HPP_INCLUDE_GUARD

#include <Optimizer/AbstractOptimizer.hpp>
#include <ObjectiveFunctional/AbstractObjectiveFunctional.hpp>

namespace controlproblem
{
    /*! \class BacktrackingOptimizer BacktrackingOptimizer.hpp
     *  \brief Class that implements the gradient method with backtracking.
     *  
     *  This class provides the implementation of the gradient method as a descent method
     *  for optimazation of funcionals. It is derived from \c AbstractOptimizer.
     *  Let \f$ J \f$ be the functional to be minimized and \f$ \psi \left( \alpha \right) \f$ the function defined as:
     *  \f[
     *      \psi \left( \alpha \right) = J \left( \mathbf{u} + \alpha\,\mathbf{d} \right)
     *  \f]
     *  with \f$ \mathbf{d} \f$ search direction.
     *  As all descent method, the generic iteration of the gradient method with backtracking can be
     *  written as:
     *  \f[
     *      \mathbf{u}_{k+1} = \mathbf{u}_k + \alpha_k\,\mathbf{d}_k\,.
     *  \f]
     *  The search direction \f$\mathbf{d}_k\f$ is chosen to be:
     *  \f[
     *      \mathbf{d}_k = - \nabla J \left(\mathbf{u}_k\right)
     *  \f]
     *  The gradient method with backtracking performs an approximated line minimization to determine the
     *  weigth \f$ \alpha \f$, imposing only the first Wolfe condition (i.e. the Armijo Rule) and
     *  using a backtracking algorithm to determine \f$ \alpha_k \f$, the weight of the 
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
    class BacktrackingOptimizer : public controlproblem::AbstractOptimizer
    {
            /************************* CONSTRUCTORS ********************/
            //! Default constructor
            BacktrackingOptimizer ();
            
            //! Constructor that overrides parameters default value
            BacktrackingOptimizer (const double& c_1, const double& alpha_0, const double& rho);

            /************************* DESTRUCTOR ********************/
            //! Destructor
            /*! Default destructor, since members of the class are trivially 
             * destructible.
             */
            virtual ~BacktrackingOptimizer () {};
            
            
            /********************** METHODS ***********************/
            //! Perform optimization on the input problem using the gradient method with backtracking
            /*! 
             *  Input arguments are:
             *  \param problem the composite differential problem that represents the primal/adjoint system
             *  \param objectiveFunctional the objective functional to be minimized
             *  \param point the starting point for the minimization algorithm
             */
            virtual void apply (const controlproblem::CompositeDifferentialProblem& problem,
                                const controlproblem::AbstractObjectiveFunctional& objectiveFunctional, 
                                dolfin::Function& point);
            
            // ---------------------------------------------------------------------------------------------//

        protected:

            // ---------------------------------------------------------------------------------------------//

        private:
            
    };
}

#endif

