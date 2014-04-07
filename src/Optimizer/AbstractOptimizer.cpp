#include <Optimizer/AbstractOptimizer.hpp>
#include <dolfin/log/dolfin_log.h>

namespace controlproblem
{
    AbstractOptimizer::AbstractOptimizer () : 
        parameters ("optimizer_parameters")
    {
        dolfin::log (dolfin::DBG, "AbstractOptimizer object created");
    }
}
