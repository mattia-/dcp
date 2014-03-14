#include <dolfin.h>
#include "Poisson.h"
#include <memory>
#include "DifferentialProblem/LinearDifferentialProblem.hpp"
#include "Utils/SubdomainType.hpp"
#include <iostream>
#include <dlfcn.h>

namespace Poisson
{
	class ExternalLoad : public dolfin::Expression
	{
		void eval (dolfin::Array<double> & values, const dolfin::Array<double> & x) const
		{
			values [0] = 10 * exp (-((x[0] - 0.5) * (x[0] - 0.5) + (x[1] - 0.5) * (x[1] - 0.5)) / 0.02);
		}
	};
	
	class NeumannCondition : public dolfin::Expression
	{
		void eval (dolfin::Array<double> & values, const dolfin::Array<double> & x) const
		{
			values [0] = sin (5 * x [0]);
		}
	};
	
	class NeumannBoundary : public dolfin::SubDomain
	{
		bool inside (const dolfin::Array<double> & x, bool on_boundary) const
		{
			return (x[1] < (0 + DOLFIN_EPS) || x[1] > (1 - DOLFIN_EPS)) && on_boundary;
		}
	};
	
	class DirichletCondition : public dolfin::Expression
	{
		void eval (dolfin::Array<double> & values, const dolfin::Array<double> & x) const
		{
			values [0] = 0;
		}
	};

	class DirichletBoundary : public dolfin::SubDomain
	{
		bool inside (const dolfin::Array<double> & x, bool on_boundary) const
		{
			return (x[0] < 0 + DOLFIN_EPS || x[0] > 1 - DOLFIN_EPS) && on_boundary;
		}
	};
	
	class UnitaryConstant : public dolfin::Expression
	{
		void eval (dolfin::Array<double>& values, const dolfin::Array<double>& x) const
		{
			values[0] = 1.0;
		}
	};
}

int main ()
{
	// declare mesh
	dolfin::UnitSquareMesh mesh (100, 100);
	
	// declare function space and linear/bilinear form
	Poisson::FunctionSpace V (mesh);
	Poisson::BilinearForm a (V, V);
	Poisson::LinearForm L (V);
	
	// declare problem specific variables
	Poisson::ExternalLoad f;
	Poisson::NeumannCondition neumannCondition;
	Poisson::DirichletCondition dirichletCondition;
	Poisson::NeumannBoundary neumannBoundary;
	Poisson::DirichletBoundary dirichletBoundary;
	
	// assign neumann condition value
	boost::shared_ptr<dolfin::Function> gg (new dolfin::Function (V));
	L.f = f;
	L.g = neumannCondition;
	Poisson::UnitaryConstant c;
	a.c = c;
	
	// set neumann condition boundary
	dolfin::FacetFunction<std::size_t> meshFacets (mesh);
	meshFacets.set_all (1);
	neumannBoundary.mark (meshFacets, 0);
	a.ds = meshFacets;
	L.ds = meshFacets;
	
	// set dirichlet condition boundary
	dolfin::DirichletBC dirichletBC (V, dirichletCondition, dirichletBoundary);

	// solve via solve() method
	dolfin::Function u (V);
	solve (a == L, u, dirichletBC);
	
	dolfin::plot (u);
	
//	for (auto i : factory.registered ())
//		std::cout << i << std::endl;
	
	// -------------------------------------------------------------------- //
	// define linear differential problem
	control_problem::LinearDifferentialProblem<Poisson::BilinearForm, Poisson::LinearForm> differentialProblem (mesh, V);

	// set solver
	differentialProblem.setSolver ("lu_solver", "default", "ilu");
	differentialProblem.setSolver ("krylov_solver", "gmres", "ilu");
	
	// set dirichlet bc
	differentialProblem.addDirichletBC (dirichletBC);
	
	boost::shared_ptr<Poisson::UnitaryConstant> c2 (new Poisson::UnitaryConstant);
	boost::shared_ptr<Poisson::ExternalLoad> f2 (new Poisson::ExternalLoad);
	boost::shared_ptr<Poisson::NeumannCondition> g2 (new Poisson::NeumannCondition);
	differentialProblem.setBilinearFormCoefficient ("c", c2);
	differentialProblem.setLinearFormCoefficient ("f", f2);
	differentialProblem.setLinearFormCoefficient ("g", g2);
	differentialProblem.setBilinearFormIntegrationSubdomains (meshFacets, control_problem::SubdomainType::BOUNDARY_FACETS);
	differentialProblem.setLinearFormIntegrationSubdomains (meshFacets, control_problem::SubdomainType::BOUNDARY_FACETS);
	
	// solve
	differentialProblem.solve ();
	differentialProblem.solve (true);
	dolfin::plot (differentialProblem.solution ());

	dolfin::interactive ();
	
	return 0;
}
