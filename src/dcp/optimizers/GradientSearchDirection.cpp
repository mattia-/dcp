#include <dcp/optimizers/GradientSearchDirection.h>

namespace dcp
{
    void GradientSearchDirection::operator() (dolfin::Function& searchDirection, const dolfin::Function gradient)
    {
        searchDirection = gradient * (-1);
    }
}
