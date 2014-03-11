#include "GenericFactory.hpp"
#include "Proxy.hpp"
#include <string>
#include <dolfin.h>

namespace control_problem
{
	typedef GenericFactory<dolfin::GenericLinearSolver, std::string> LinearSolverFactory;
}

