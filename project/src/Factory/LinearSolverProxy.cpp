#include <Factory/LinearSolverFactory.hpp>

// this file will be compiled into a library linked to the library generated from LinearDifferentialSolver.
// Once this library is loaded, the objects defined where will be created, and the corresponding creation rules will
// be added to the factory
namespace
{
	control_problem::Proxy <control_problem::LinearSolverFactory, dolfin::LUSolver> LU_S ("lu_solver");
	control_problem::Proxy <control_problem::LinearSolverFactory, dolfin::KrylovSolver> K_S ("krylov_solver");
}
