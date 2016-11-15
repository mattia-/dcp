#include <dcp/optimizers/GradientSearchDirection.h>

namespace dcp
{
    dolfin::Function GradientSearchDirection::operator() (const dolfin::Function& gradient)
    {
        // use temporary object because dolfin does not provide a constructor for dolfin::Function based on
        // dolfin::FunctionAXPY
        dolfin::Function result (gradient);
        result = gradient * (-1);
        return result;
    }



    dcp::TimeDependentFunction GradientSearchDirection::operator() (const dcp::TimeDependentFunction& gradient)
    {
        return gradient * (-1);
    }
}
