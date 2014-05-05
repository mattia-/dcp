#ifndef SRC_OPTIMIZER_DIRICHLETCONTROLVALUEUPDATER_HPP_INCLUDE_GUARD
#define SRC_OPTIMIZER_DIRICHLETCONTROLVALUEUPDATER_HPP_INCLUDE_GUARD

#include <DifferentialProblem/CompositeDifferentialProblem.hpp>
#include <dolfin/function/GenericFunction.h>
#include <dolfin/function/FunctionSpace.h>
#include <dolfin/mesh/SubDomain.h>
#include <string>
#include <boost/shared_ptr.hpp>

namespace DCP
{
    /*! \class DirichletControlValueUpdater DirichletControlValueUpdater.hpp
     *  \brief Class to update the value of the control variable in Dirichlet boundary control problems.
     *  
     *  This class is a functor which can be passed to the method \c apply() of any class
     *  of the \c AbstractOptimizer hierarchy, which will use it to update the value of
     *  the control parameter in the \c DifferentialProblem (also passed to the method \c apply()
     *  of the same class) as the optimization proceeds.
     */
    class DirichletControlValueUpdater
    {
        // ---------------------------------------------------------------------------------------------//

        public:
            /************************* CONSTRUCTORS ********************/
            //! Default constructor is deleted
            DirichletControlValueUpdater () = delete;
            
            //! Constructor
            /*! 
             *  Input arguments are:
             *  \param problemName string that identifies the problem (in the \c CompositeDifferentialProblem object 
             *  passed as input to <tt> this->operator() ()</tt> ) which contains the control parameter to be updated
             *  \param dirichletBCName the name under which the control Dirichlet boundary condition is stored in the 
             *  protected member \c map in the problem
             *  \param dirichletBoundary the boundary over which the control Dirichlet condition should be enforced
             *  \param functionSpace the function space (possibly a subspace) on which the dirichlet condition should be enforced.
             */
            DirichletControlValueUpdater (const std::string& problemName, 
                                          const std::string& dirichletBCName,
                                          const dolfin::SubDomain& dirichletBoundary,
                                          boost::shared_ptr<const dolfin::FunctionSpace> functionSpace);


            /************************* OPERATORS ********************/
            //! Call operator. This will actually perform the updating of the control parameter
            /*! 
             *  Input parametes are:
             *  \param compositeProblem the problem on which to operate
             *  \param dirichletBCValue the new value for the control Dirichlet boundary condition 
             */
            void operator() (DCP::CompositeDifferentialProblem& compositeProblem, 
                             const dolfin::GenericFunction& dirichletBCValue) const;

            // ---------------------------------------------------------------------------------------------//

        protected:
            //! The name identifying the problem inside the \c CompositeDifferentialProblem which contains the
            //! control parameter
            std::string problemName_;

            //! The name under which the control Dirichlet boundary condition is stored in the 
            //! protected member \c map in the problem
            std::string dirichletBCName_;

            //! The boundary over which the control Dirichlet condition should be enforced
            const dolfin::SubDomain& dirichletBoundary_;

            //! The function space (possibly a subspace) on which the dirichlet condition should be enforced.
            boost::shared_ptr<const dolfin::FunctionSpace> functionSpace_;
    };
}

#endif


