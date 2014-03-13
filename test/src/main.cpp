#include <dolfin.h>
#include "Poisson.h"
#include <memory>
#include "DifferentialProblem/LinearDifferentialProblem.hpp"

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
		
		public:
		~DirichletBoundary ()
		{
			std::cout << "I am destroyed :)" << std::endl;
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
	
	dolfin::info (mesh, false);
	// declare function space and linear/bilinear form
	Poisson::FunctionSpace V (mesh);
	Poisson::BilinearForm a (V, V);
	Poisson::LinearForm L (V);
	dolfin::cout << "ciao" << dolfin::endl;
	std::unique_ptr<Poisson::BilinearForm> foo (new Poisson::BilinearForm (a));
	dolfin::cout << "ciao" << dolfin::endl;
//	foo.reset (&a);
	dolfin::cout << "ciao" << dolfin::endl;
//	foo.reset ();
	if (foo == nullptr)
		dolfin::cout << "ciao NULL" << dolfin::endl;
	
	// declare problem specific variables
	Poisson::ExternalLoad f;
	Poisson::NeumannCondition neumannCondition;
	Poisson::DirichletCondition dirichletCondition;
	Poisson::NeumannBoundary neumannBoundary;
	Poisson::DirichletBoundary dirichletBoundary;
	
	// assign neumann condition value
	boost::shared_ptr<dolfin::Function> gg (new dolfin::Function (V));
	L.set_coefficient ("g", gg);
	L.f = f;
//	L.g = neumannCondition;
	boost::shared_ptr<Poisson::UnitaryConstant> c (new Poisson::UnitaryConstant);
//	Poisson::UnitaryConstant c;
//	a.c = *c;
	a.set_coefficient ("c", c);
	decltype(a) ba (a);
	typedef decltype(a) myType;
	ba.set_coefficient ("c", c);
	foo->set_coefficient(0, c);
	
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
	
	// solve with LU solver
	dolfin::Function u2 (V);
	boost::shared_ptr<dolfin::Matrix> A2 (new dolfin::Matrix);
	dolfin::Vector b2;
	
	dolfin::assemble (*A2, a);
	dolfin::assemble (b2, L);
	
	dirichletBC.apply (*A2, b2);
	
	dolfin::LUSolver solver;
	std::unique_ptr<dolfin::GenericLinearSolver> solver2 (new dolfin::LUSolver (solver));
	solver2->set_operator (A2);
	solver2->solve (*u2.vector (), b2);
	dolfin::plot (u2);
	
	// solve with Krylov solver
	dolfin::Function u3 (V);
	dolfin::Matrix A3;
//	boost::shared_ptr<dolfin::Matrix> A3 (new dolfin::Matrix);
	dolfin::Vector b3;
	*gg = neumannCondition;	
	dolfin::assemble (A3, a);
//	dolfin::assemble (*A3, a);
	dolfin::assemble (b3, L);
	
	dirichletBC.apply (A3, b3);
//	dirichletBC.apply (*A3, b3);
	
	dolfin::KrylovSolver kSolver ("gmres", "ilu");
	kSolver.set_operator (boost::shared_ptr<dolfin::Matrix> (new dolfin::Matrix (A3)));
	kSolver.solve (*u3.vector (), b3);
	
	dolfin::plot (u3);
	
	dolfin::interactive ();
	
	std::vector <int> vec (10, 5);
	for (auto i : vec)
	{
		std::cout << i << std::endl;
	}
		
	std::unique_ptr <int> intero (new int);
	*intero = 8;
	std::cout << *intero << std::endl;
	
	{
	std::unique_ptr <Poisson::DirichletBoundary> aa (new Poisson::DirichletBoundary);
	}
	
	int n;
	std::cin >> n;

	return 0;
}
