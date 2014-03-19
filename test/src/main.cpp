#include "Poisson.h"
#include "NavierStokes.h"
#include <memory>
#include "DifferentialProblem/LinearDifferentialProblem.hpp"
#include "DifferentialProblem/NonLinearDifferentialProblem.hpp"
#include "Utils/SubdomainType.hpp"
#include <iostream>
#include <dlfcn.h>
#include <dolfin.h>

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

namespace NavierStokes
{
	class MovingLidBoundary : public dolfin::SubDomain
	{
		bool inside (const dolfin::Array<double> & x, bool on_boundary) const
		{
			return x [1] > (1.0 - DOLFIN_EPS) && on_boundary;
		}
	};

	class NoSlipBoundary : public dolfin::SubDomain
	{
		bool inside (const dolfin::Array<double> & x, bool on_boundary) const
		{
			return x [1] < (1.0 - DOLFIN_EPS) && on_boundary;
		}
	};
}

int main ()
{
	// ************************************************************************************************ //
	// *********************************** LINEAR PROBLEM TESTS *************************************** //
	// ************************************************************************************************ //
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
//	std::shared_ptr<dolfin::UnitSquareMesh> mesh2 (new dolfin::UnitSquareMesh (100, 100));
//	std::shared_ptr<Poisson::FunctionSpace> fs2 (new Poisson::FunctionSpace (*mesh2));
//	control_problem::LinearDifferentialProblem<Poisson::BilinearForm, Poisson::LinearForm> differentialProblem (mesh, V);
	control_problem::LinearDifferentialProblem<Poisson::BilinearForm, Poisson::LinearForm> differentialProblem
		(std::move (mesh), std::move (V));

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
	
	// -------------------------------------------------------------------- //
	// define linear differential problem
//	std::shared_ptr<dolfin::UnitSquareMesh> mesh2 (new dolfin::UnitSquareMesh (100, 100));
//	std::shared_ptr<Poisson::FunctionSpace> fs2 (new Poisson::FunctionSpace (*mesh2));
	dolfin::UnitSquareMesh mesh2 (100, 100);
	Poisson::FunctionSpace fs2 (mesh2);
	control_problem::LinearDifferentialProblem<Poisson::BilinearForm, Poisson::LinearForm> differentialProblem2 (mesh2, fs2);

	// set solver
	differentialProblem2.setSolver ("lu_solver", "default", "ilu");
	differentialProblem2.setSolver ("krylov_solver", "gmres", "ilu");
	
	// set dirichlet bc
//	dolfin::DirichletBC dirichletBC2 (*fs2, dirichletCondition, dirichletBoundary);
	dolfin::DirichletBC dirichletBC2 (fs2, dirichletCondition, dirichletBoundary);
	differentialProblem2.addDirichletBC (dirichletBC2);
	
	boost::shared_ptr<Poisson::UnitaryConstant> c3 (new Poisson::UnitaryConstant);
	boost::shared_ptr<Poisson::ExternalLoad> f3 (new Poisson::ExternalLoad);
	boost::shared_ptr<Poisson::NeumannCondition> g3 (new Poisson::NeumannCondition);
	differentialProblem2.setBilinearFormCoefficient ("c", c3);
	differentialProblem2.setLinearFormCoefficient ("f", f3);
	differentialProblem2.setLinearFormCoefficient ("g", g3);
	

	dolfin::FacetFunction<std::size_t> meshFacets2 (mesh2);
	meshFacets2.set_all (1);
	neumannBoundary.mark (meshFacets2, 0);
	differentialProblem2.setBilinearFormIntegrationSubdomains (meshFacets2, control_problem::SubdomainType::BOUNDARY_FACETS);
	differentialProblem2.setLinearFormIntegrationSubdomains (meshFacets2, control_problem::SubdomainType::BOUNDARY_FACETS);
	
	// solve
	differentialProblem2.solve ();
	differentialProblem2.solve (true);
	dolfin::plot (differentialProblem2.solution ());

//	dolfin::interactive ();
	
	// -----------------------------------------------
	
	// ************************************************************************************************ //
	// ********************************* NON LINEAR PROBLEM TESTS ************************************* //
	// ************************************************************************************************ //
	
	// define mesh
	dolfin::UnitSquareMesh NLmesh (20, 20);
	
	// define vector space
	NavierStokes::FunctionSpace NLV (NLmesh);
	
	// define variational form
	NavierStokes::ResidualForm NLF (NLV);

	// define coefficients
	dolfin::Constant NLnu (1e-6);
	NLF.nu = NLnu;
	
	// define NLsolution
	dolfin::Function NLsolution (NLV);
	
	// boundary condition
	std::vector<const dolfin::DirichletBC*> NLboundaryConditions;
	
	// no-slip boundary conditions
	NavierStokes::NoSlipBoundary NLnoSlipBoundary;
	dolfin::Constant NLnoSlipValue (0.0, 0.0);
	dolfin::DirichletBC NLnoSlipCondition (*NLV[0], NLnoSlipValue, NLnoSlipBoundary);
	NLboundaryConditions.emplace_back (&NLnoSlipCondition);
	
	// moving lid boundary conditions
	dolfin::DirichletBC NLmovingLidCondition (NLV[0], boost::shared_ptr<dolfin::Constant> (new dolfin::Constant (1.0, 0.0)), 
	                                        boost::shared_ptr<NavierStokes::MovingLidBoundary> (new NavierStokes::MovingLidBoundary));
	NLboundaryConditions.emplace_back (&NLmovingLidCondition);
	boost::shared_ptr<dolfin::Function> NLinitialGuess (new dolfin::Function (NLsolution));
	NLF.set_coefficient ("trial", NLinitialGuess);
	
	NavierStokes::JacobianForm NLJ (NLV, NLV);
	{
	boost::shared_ptr<dolfin::Function> NLinitialGuess2 (new dolfin::Function (NLsolution));
	NLJ.set_coefficient ("trial", NLinitialGuess2);
	}
	boost::shared_ptr<dolfin::Constant> NLnup (new dolfin::Constant (1e-6));
	NLJ.set_coefficient ("nu", NLnup);	
	
	solve (NLF == 0, NLsolution, NLboundaryConditions, NLJ);
	
	dolfin::plot (NLsolution[0]);
	dolfin::plot (NLsolution[0][0]);
	dolfin::plot (NLsolution[0][1]);
	dolfin::plot (NLsolution[1]);
	
	// --------------------------------------------------------

	boost::shared_ptr<dolfin::Function> NLinitialGuess3 (new dolfin::Function (NLV));
	NLF.set_coefficient ("trial", NLinitialGuess3);
	
	{
	boost::shared_ptr<dolfin::Function> NLinitialGuess2 (new dolfin::Function (NLV));
	NLJ.set_coefficient ("trial", NLinitialGuess2);
	}
	
	solve (NLF == 0, NLsolution, NLboundaryConditions, NLJ);
	
	dolfin::plot (NLsolution[0]);
	dolfin::plot (NLsolution[0][0]);
	dolfin::plot (NLsolution[0][1]);
	dolfin::plot (NLsolution[1]);
	
	// --------------------------------------------------------
	
	dolfin::begin ("Now with NonLinearDifferentialProblem class");
		dolfin::begin ("First solution");
		
			// define mesh
			dolfin::UnitSquareMesh NLmesh2 (20, 20);
			
			// define vector space
			NavierStokes::FunctionSpace NLV2 (NLmesh2);
			
			control_problem::NonLinearDifferentialProblem<NavierStokes::ResidualForm, NavierStokes::JacobianForm> 
				NLdifferentialProblem (NLmesh2, NLV2, "trial");

			// define coefficients
			boost::shared_ptr<dolfin::Constant> NLnu2 (new dolfin::Constant (1e-6));
			NLdifferentialProblem.setResidualFormCoefficient ("nu", NLnu2);
			
			
			// no-slip boundary conditions
			NavierStokes::NoSlipBoundary NLnoSlipBoundary2;
			dolfin::Constant NLnoSlipValue2 (0.0, 0.0);
			dolfin::DirichletBC NLnoSlipCondition2 (*NLV2[0], NLnoSlipValue2, NLnoSlipBoundary2);
			NLdifferentialProblem.addDirichletBC (NLnoSlipCondition2);
			
			// moving lid boundary conditions
			dolfin::DirichletBC NLmovingLidCondition2 
				(NLV[0], boost::shared_ptr<dolfin::Constant> (new dolfin::Constant (1.0, 0.0)), 
				 boost::shared_ptr<NavierStokes::MovingLidBoundary> (new NavierStokes::MovingLidBoundary));
			NLdifferentialProblem.addDirichletBC (NLmovingLidCondition2);
			
		//	boost::shared_ptr<dolfin::Function> NLinitialGuess2 (new dolfin::Function (NLsolution));
			
			NLdifferentialProblem.setJacobianFormCoefficient ("nu", NLnu2);	
			
			NLdifferentialProblem.solve ();
		
		dolfin::end ();
		
		dolfin::begin ("Setting initial guess and parameters");
		
			dolfin::Function init (NLV);
		
//			NLdifferentialProblem.setInitialGuess (NLdifferentialProblem.solution ());
			NLdifferentialProblem.setInitialGuess (init);

			dolfin::Parameters p ("my_parameters");
			dolfin::Parameters q ("newton_solver");
			q.add ("linear_solver", "lu");
			p.add (q);
			
			NLdifferentialProblem.setSolverParameters (p);	
			dolfin::info (NLdifferentialProblem.solverParameters (), true);
			NLdifferentialProblem.solve (p);
			
			dolfin::end ();
	
	dolfin::end ();
	
	dolfin::plot (NLdifferentialProblem.solution ()[0]);
	dolfin::plot (NLdifferentialProblem.solution ()[0][0]);
	dolfin::plot (NLdifferentialProblem.solution ()[0][1]);
	dolfin::plot (NLdifferentialProblem.solution ()[1]);

	dolfin::interactive ();

	return 0;
}
