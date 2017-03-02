#include <dcp/utils/plots.h>

namespace dolfin
{
    void plot (dcp::TimeDependentFunction& function, std::string title, std::string mode)
    {
        if (mode == "pause")
        {
            function.plot (title, true);
        }
        else if (mode == "continue")
        {
            function.plot (title, false);
        }
        else
        {
            dolfin::dolfin_error ("dcp: plots.cpp",
                                  "plot object of type dcp::TimeDependentFunction",
                                  "Unknown plot mode \"%s\"",
                                  mode.c_str ());
        }
    }



    void plot (const std::shared_ptr<dcp::TimeDependentFunction> function, std::string title, std::string mode)
    {
        if (mode == "pause")
        {
            function->plot (title, true);
        }
        else if (mode == "continue")
        {
            function->plot (title, false);
        }
        else
        {
            dolfin::dolfin_error ("dcp: plots.cpp",
                                  "plot object of type dcp::TimeDependentFunction",
                                  "Unknown plot mode \"%s\"",
                                  mode.c_str ());
        }
    }
}
