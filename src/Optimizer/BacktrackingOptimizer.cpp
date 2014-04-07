#include <Optimizer/BacktrackingOptimizer.hpp>
#include <dolfin/parameter/Parameters.h>
#include <dolfin/log/dolfin_log.h>

namespace controlproblem
{
    BacktrackingOptimizer::BacktrackingOptimizer () : 
        AbstractOptimizer ()
    {
        dolfin::begin (dolfin::DBG, "Creating BacktrackingOptimizer object...");
        
        dolfin::log (dolfin::DBG, "Setting up parameters...");
        parameters.add ("descent_method", "backtracking_gradient_method");
        parameters.add ("c_1", 1e-3);
        parameters.add ("alpha_0", 0.5);
        parameters.add ("rho", 0.5);
        
        dolfin::log (dolfin::DBG, "BacktrackingOptimizer object created");
        
        dolfin::end ();
    }
    
    
    BacktrackingOptimizer::BacktrackingOptimizer (const double& c_1, const double& alpha_0, const double& rho) :
        AbstractOptimizer ()
    {
        dolfin::begin (dolfin::DBG, "Creating BacktrackingOptimizer object...");
        
        dolfin::log (dolfin::DBG, "Setting up parameters...");
        parameters.add ("descent_method", "backtracking_gradient_method");
        parameters.add ("c_1", c_1);
        parameters.add ("alpha_0", alpha_0);
        parameters.add ("rho", rho);
        
        dolfin::log (dolfin::DBG, "BacktrackingOptimizer object created");
        
        dolfin::end ();
    }
    
    


    void BacktrackingOptimizer::apply (const controlproblem::CompositeDifferentialProblem& problem,
                                       const controlproblem::AbstractObjectiveFunctional& objectiveFunctional, 
                                       dolfin::Function& point)
    {
       // get parameters values
       double c_1     = this->parameters ["c_1"];
       double alpha_0 = this->parameters ["alpha_0"];
       double rho     = this->parameters ["rho"];
       
       double alpha = alpha_0;
       
       // lambra functions to check if Armijo condition is satisfied
       auto psi = [&] () {return objectiveFunctional.evaluateFunctional ()<++>}<++>
       auto armijoConditionCheck 

        
    }



}

