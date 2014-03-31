#include "Poisson.h"
#include "NavierStokes.h"
#include <memory>
#include "DifferentialProblem/LinearDifferentialProblem.hpp"
#include "DifferentialProblem/NonlinearDifferentialProblem.hpp"
#include "DifferentialProblem/CompositeDifferentialProblem.hpp"
#include "Utils/SubdomainType.hpp"
#include <iostream>
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
    // =============================================================================================================== //
    // =============================================================================================================== //
    // =============================================================================================================== //
    
    dolfin::set_log_level (dolfin::DBG);
//    dolfin::set_log_level (dolfin::PROGRESS);
    
    control_problem::CompositeDifferentialProblem cdp;
    
    // declare meshes
    std::shared_ptr<dolfin::UnitSquareMesh> mesh (new dolfin::UnitSquareMesh (100, 100));
    std::shared_ptr<dolfin::UnitSquareMesh> NLmesh (new dolfin::UnitSquareMesh (20, 20));
    
    // declare function spaces
    std::shared_ptr<Poisson::FunctionSpace> V (new Poisson::FunctionSpace (*mesh));
    std::shared_ptr<NavierStokes::FunctionSpace> NLV (new NavierStokes::FunctionSpace (*NLmesh));
    
    // declare problems
    control_problem::LinearDifferentialProblem<Poisson::BilinearForm, Poisson::LinearForm> ldp (*mesh, *V, "lu_solver"); 
//    std::unique_ptr<control_problem::AbstractDifferentialProblem> ldp 
//        (new control_problem::LinearDifferentialProblem<Poisson::BilinearForm, Poisson::LinearForm> (*mesh, *V, "lu_solver"));
    control_problem::NonlinearDifferentialProblem<NavierStokes::ResidualForm, NavierStokes::JacobianForm> nldp (*NLmesh, *NLV, "trial");
//    std::unique_ptr<control_problem::AbstractDifferentialProblem> nldp 
//        (new control_problem::NonlinearDifferentialProblem<NavierStokes::ResidualForm, NavierStokes::JacobianForm> (*NLmesh, *NLV, "trial"));
    
    cdp.addProblem ("ldp", ldp);
    cdp.addProblem ("nldp", nldp);
//    cdp.print ();
  
    std::vector<std::string> v;
    v.push_back ("ldp");
    v.push_back ("nldp");
    v.push_back ("ldp");
    v.push_back ("nldp");
    cdp.reorderProblems (v);
//    cdp.print ();
    
    v.pop_back ();
    cdp.reorderProblems (v);
//    cdp.print ();
  
  // -----------------------------------------------------------------------------------------------------//
  // --------------------------------------- LINEAR PROBLEM ----------------------------------------------//
  // -----------------------------------------------------------------------------------------------------//
    
    Poisson::NeumannCondition neumannCondition;
    Poisson::DirichletCondition dirichletCondition;
    Poisson::NeumannBoundary neumannBoundary;
    Poisson::DirichletBoundary dirichletBoundary;
    
    boost::shared_ptr<dolfin::FacetFunction<std::size_t>> meshFacets (new dolfin::FacetFunction<std::size_t> (*mesh));
    meshFacets->set_all (1);
    neumannBoundary.mark (*meshFacets, 0);
    dolfin::DirichletBC dirichletBC (*V, dirichletCondition, dirichletBoundary);

    
//    boost::shared_ptr<Poisson::UnitaryConstant> c2 (new Poisson::UnitaryConstant);
    Poisson::UnitaryConstant c2;
    boost::shared_ptr<Poisson::ExternalLoad> f2 (new Poisson::ExternalLoad);
    boost::shared_ptr<Poisson::NeumannCondition> g2 (new Poisson::NeumannCondition);
    cdp["ldp"].setCoefficient ("bilinear_form", dolfin::reference_to_no_delete_pointer (c2), "c");
    cdp["ldp"].setCoefficient ("linear_form", f2, "f");
    cdp["ldp"].setCoefficient ("linear_form", g2, "g");
    
    cdp["ldp"].setIntegrationSubdomains ("bilinear_form", meshFacets, control_problem::SubdomainType::BOUNDARY_FACETS);
    cdp["ldp"].setIntegrationSubdomains ("linear_form", meshFacets, control_problem::SubdomainType::BOUNDARY_FACETS);
    
    cdp["ldp"].addDirichletBC (dirichletBC);
    
    // solve
    cdp.solve ("ldp");
//    dolfin::plot (cdp["ldp"].solution ());

    // -----------------------------------------------------------------------------------------------------//
    // ------------------------------------- NON LINEAR PROBLEM --------------------------------------------//
    // -----------------------------------------------------------------------------------------------------//
     
    // boundary condition
    std::vector<const dolfin::DirichletBC*> NLboundaryConditions;
    
    // no-slip boundary conditions
    NavierStokes::NoSlipBoundary NLnoSlipBoundary;
    dolfin::Constant NLnoSlipValue (0.0, 0.0);
    dolfin::DirichletBC NLnoSlipCondition (*((*NLV)[0]), NLnoSlipValue, NLnoSlipBoundary);
    NLboundaryConditions.emplace_back (&NLnoSlipCondition);
    
    // moving lid boundary conditions
    dolfin::DirichletBC NLmovingLidCondition ((*NLV)[0], boost::shared_ptr<dolfin::Constant> (new dolfin::Constant (1.0, 0.0)), 
                                            boost::shared_ptr<NavierStokes::MovingLidBoundary> (new NavierStokes::MovingLidBoundary));
    NLboundaryConditions.emplace_back (&NLmovingLidCondition);
    
//    boost::shared_ptr<dolfin::Function> NLinitialGuess (new dolfin::Function (NLsolution));
    
    boost::shared_ptr<dolfin::Constant> NLnu2 (new dolfin::Constant (1e-6));
    cdp["nldp"].setCoefficient ("residual_form", NLnu2, "nu");
    cdp["nldp"].setCoefficient ("jacobian_form", NLnu2, "nu");
    
    for (auto i : NLboundaryConditions)
    {
        cdp["nldp"].addDirichletBC (*i);
    }
    
    cdp.solve ("nldp", true);
//    dolfin::plot (cdp.problem ("nldp").solution ()[0]);
//    dolfin::plot (cdp.problem ("nldp").solution ()[0][0]);
//    dolfin::plot (cdp.problem ("nldp").solution ()[0][1]);
//    dolfin::plot (cdp.problem ("nldp").solution ()[1]);
    
//    dolfin::cout << dolfin::endl;
//    dolfin::cout << dolfin::endl;
//    dolfin::cout << dolfin::endl;
//    dolfin::cout << dolfin::endl;
//    dolfin::cout << "-----------------------------------------------------------------" << dolfin::endl;
//    dolfin::cout << dolfin::endl;
//    dolfin::cout << dolfin::endl;
//    dolfin::cout << dolfin::endl;
//    dolfin::cout << dolfin::endl;
    
    cdp.solve ();
    cdp.solve (true);
    
    dolfin::plot (cdp.solution ("nldp")[0]);
    dolfin::plot (cdp.solution ("nldp")[0][0]);
    dolfin::plot (cdp.solution ("nldp")[1]);
    dolfin::plot (cdp[1].solution ()[0][1]);
    dolfin::plot (cdp.solution ("ldp"));
    dolfin::interactive ();
    
    dolfin::cout << "Problem size is: " << cdp.size () << dolfin::endl;
    return 0;
}
