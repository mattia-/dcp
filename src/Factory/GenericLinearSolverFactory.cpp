#include "GenericFactory.hpp"
#include "Proxy.hpp"
#include <string>
#include <dolfin.h>

// this file will be compiled into a library linked to the library generated from LinearDifferentialSolver.
// Once this library is loaded, the objects defined where will be created, and the corresponding creation rules will
// be added to the factory
namespace
{
	Proxy <GenericFactory <dolfin::GenericLinearSolver, std::string>, dolfin::LUSolver> LU_S ("lu_solver");
	Proxy <GenericFactory <dolfin::GenericLinearSolver, std::string>, dolfin::krylovSolver> K_S ("krylov_solver");
}
