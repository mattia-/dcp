#ifndef HH__LINEARSOLVERFACTORY__HH
#define HH__LINEARSOLVERFACTORY__HH
#include <Factory/GenericFactory.hpp>
#include <Factory/Proxy.hpp>
#include <string>
#include <dolfin.h>

namespace control_problem
{
	typedef GenericFactory <dolfin::GenericLinearSolver, std::string> LinearSolverFactory;
}
#endif
