#include <dcp/utils/plots.h>

namespace dolfin
{
    void plot (dcp::TimeDependentFunction& function, std::string title, const bool& pause)
    {
        function.plot (title, pause);
    }



    void plot (const std::shared_ptr<dcp::TimeDependentFunction> function, std::string title, const bool& pause)
    {
        function->plot (title, pause);
    }
}
