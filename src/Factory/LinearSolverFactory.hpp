#ifndef SRC_FACTORY_LINEARSOLVERFACTORY_HPP_INCLUDE_GUARD
#define SRC_FACTORY_LINEARSOLVERFACTORY_HPP_INCLUDE_GUARD
#include <Factory/GenericFactory.hpp>
#include <Factory/Proxy.hpp>
#include <string>
#include <dolfin.h>

namespace controlproblem
{
    typedef GenericFactory <dolfin::GenericLinearSolver, std::string> LinearSolverFactory;
}
#endif
