#include <Optimizer/DirichletControlValueUpdater.hpp>
#include <DifferentialProblem/AbstractDifferentialProblem.hpp>
#include <dolfin/fem/DirichletBC.h>

namespace controlproblem
{
    /************************* CONSTRUCTORS ********************/
    DirichletControlValueUpdater::DirichletControlValueUpdater (const std::string& problemName, 
                                                                const std::string& dirichletBCName,
                                                                const dolfin::SubDomain& dirichletBoundary,
                                                                boost::shared_ptr<const dolfin::FunctionSpace> functionSpace) : 
        problemName_ (problemName),
        dirichletBCName_ (dirichletBCName),
        dirichletBoundary_ (dirichletBoundary),
        functionSpace_ (functionSpace)
    {  }



    /************************* OPERATORS ********************/
    void DirichletControlValueUpdater::operator() (controlproblem::CompositeDifferentialProblem& compositeProblem, 
                                                   const dolfin::GenericFunction& dirichletBCValue) const
    {
        controlproblem::AbstractDifferentialProblem& problem = compositeProblem [problemName_];
        
        problem.removeDirichletBC (dirichletBCName_);
        
        problem.addDirichletBC (dolfin::DirichletBC (*(functionSpace_), dirichletBCValue, dirichletBoundary_),
                                dirichletBCName_);
    }
}
