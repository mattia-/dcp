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
    
    dolfin::cout << "************************************************************************************************" << dolfin::endl;
    dolfin::cout << "*********************************** LINEAR PROBLEM TESTS ***************************************" << dolfin::endl;
    dolfin::cout << "************************************************************************************************" << dolfin::endl;
    dolfin::cout << "------------------------------------------------------------------------------------------------" << dolfin::endl;
    dolfin::cout << "----------------------------------- NORMAL SOLUTION PROCESS ------------------------------------" << dolfin::endl;
    dolfin::cout << "------------------------------------------------------------------------------------------------" << dolfin::endl;
    
    // declare mesh
    boost::shared_ptr<dolfin::UnitSquareMesh> mesh (new dolfin::UnitSquareMesh (100, 100));
    
    // declare function space and linear/bilinear form
    boost::shared_ptr<Poisson::FunctionSpace> V (new Poisson::FunctionSpace (*mesh));
    Poisson::BilinearForm a (*V, *V);
    Poisson::LinearForm L (*V);
    
    // declare problem specific variables
    Poisson::ExternalLoad f;
    Poisson::NeumannCondition neumannCondition;
    Poisson::DirichletCondition dirichletCondition;
    Poisson::NeumannBoundary neumannBoundary;
    Poisson::DirichletBoundary dirichletBoundary;
    
    // assign neumann condition value
    boost::shared_ptr<dolfin::Function> gg (new dolfin::Function (*V));
    L.f = f;
    L.g = neumannCondition;
    Poisson::UnitaryConstant c;
    a.c = c;
    
    // set neumann condition boundary
    dolfin::FacetFunction<std::size_t> meshFacets (*mesh);
    meshFacets.set_all (1);
    neumannBoundary.mark (meshFacets, 0);
    a.ds = meshFacets;
    L.ds = meshFacets;
    
    // set dirichlet condition boundary
    dolfin::DirichletBC dirichletBC (*V, dirichletCondition, dirichletBoundary);

    // solve via solve() method
    dolfin::Function u (*V);
    solve (a == L, u, dirichletBC);
    
    dolfin::plot (u);
    
    dolfin::set_log_level (dolfin::DBG);
    
    dolfin::cout << "------------------------------------------------------------------------------------------------" << dolfin::endl;
    dolfin::cout << "------------------------------ WITH LINEARDIFFERENTIALPROBLEM CLASS ----------------------------" << dolfin::endl;
    dolfin::cout << "------------------------------------------------------------------------------------------------" << dolfin::endl;
    dolfin::cout << dolfin::endl;
    dolfin::cout << "....................... FIRST CONSTRUCTOR ......................." << dolfin::endl;
    
    // define linear differential problem
    DCP::LinearDifferentialProblem<Poisson::BilinearForm, Poisson::LinearForm> differentialProblem (*mesh, *V, "lu_solver");

    // set solver
//    dolfin::info (differentialProblem.parameters, true);
    differentialProblem.parameters ["desired_solver_type"] = "krylov_solver";
    differentialProblem.parameters ["desired_solver_method"] = "gmres";
    differentialProblem.parameters ["desired_solver_preconditioner"] = "ilu";
//    dolfin::info (differentialProblem.parameters, true);
    
    differentialProblem.update ();
//    dolfin::info (differentialProblem.parameters, true);

    // set dirichlet bc
    differentialProblem.addDirichletBC (dirichletBC, "my_name");
    
    boost::shared_ptr<Poisson::UnitaryConstant> c2 (new Poisson::UnitaryConstant);
    boost::shared_ptr<Poisson::ExternalLoad> f2 (new Poisson::ExternalLoad);
    boost::shared_ptr<Poisson::NeumannCondition> g2 (new Poisson::NeumannCondition);
    differentialProblem.setCoefficient ("bilinear_form", c2, "c");
    differentialProblem.setCoefficient ("linear_form", f2, "f");
    differentialProblem.setCoefficient ("linear_form", g2, "g");
    
    differentialProblem.setIntegrationSubdomains ("bilinear_form",
                                                  dolfin::reference_to_no_delete_pointer (meshFacets),
                                                  DCP::SubdomainType::BOUNDARY_FACETS);
    differentialProblem.setIntegrationSubdomains ("linear_form", 
                                                  dolfin::reference_to_no_delete_pointer (meshFacets),
                                                  DCP::SubdomainType::BOUNDARY_FACETS);
    
    // solve
    differentialProblem.solve ();
    differentialProblem.solve (true);
    differentialProblem.solve ();
    differentialProblem.parameters ["force_reassemble_system"] = true;
    differentialProblem.solve ();
    differentialProblem.parameters ["force_reassemble_system"] = false;
    differentialProblem.solve ();
    differentialProblem.parameters ["desired_solver_type"] = "lu_solver";
    differentialProblem.parameters ["desired_solver_method"] = "default";
    differentialProblem.solve ();
    dolfin::plot (differentialProblem.solution ());
    
    dolfin::Form a1 (differentialProblem.bilinearForm ());
    dolfin::Form a2 (differentialProblem.linearForm ());
    dolfin::Matrix a4 (differentialProblem.linearOperator ());
    dolfin::Vector a5 (differentialProblem.rhs ());
    
    // ............................................................................................................ //
    dolfin::cout << dolfin::endl;
    dolfin::cout << "....................... SECOND CONSTRUCTOR ......................." << dolfin::endl;
    // define linear differential problem
    DCP::LinearDifferentialProblem<Poisson::BilinearForm, Poisson::LinearForm> differentialProblem2 (mesh, V);

//    dolfin::info (differentialProblem2.parameters, true);
    differentialProblem2.parameters ["desired_solver_type"] = "krylov_solver";
    differentialProblem2.parameters ["desired_solver_method"] = "gmres";
    differentialProblem2.parameters ["desired_solver_preconditioner"] = "ilu";
//    dolfin::info (differentialProblem2.parameters, true);
    
    differentialProblem2.update ();
//    dolfin::info (differentialProblem2.parameters, true);

    // set dirichlet bc
    differentialProblem2.addDirichletBC (dirichletBC);
    
    differentialProblem2.setCoefficient ("bilinear_form", c2, "c");
    differentialProblem2.setCoefficient ("linear_form", f2, "f");
    differentialProblem2.setCoefficient ("linear_form", g2, "g");
    
    differentialProblem2.setIntegrationSubdomains ("bilinear_form", 
                                                   dolfin::reference_to_no_delete_pointer (meshFacets),
                                                   DCP::SubdomainType::BOUNDARY_FACETS);
    differentialProblem2.setIntegrationSubdomains ("linear_form", 
                                                   dolfin::reference_to_no_delete_pointer (meshFacets),
                                                   DCP::SubdomainType::BOUNDARY_FACETS);
    
    // solve
    differentialProblem2.solve ();
    differentialProblem2.solve (true);
    differentialProblem2.solve ();
    differentialProblem2.parameters ["force_reassemble_system"] = true;
    differentialProblem2.solve ();
    differentialProblem2.parameters ["force_reassemble_system"] = false;
    differentialProblem2.solve ();
    differentialProblem2.parameters ["desired_solver_type"] = "lu_solver";
    differentialProblem2.parameters ["desired_solver_method"] = "default";
    differentialProblem2.solve ();
    dolfin::plot (differentialProblem2.solution ());
    
    dolfin::Form a11 (differentialProblem2.bilinearForm ());
    dolfin::Form a21 (differentialProblem2.linearForm ());
    dolfin::Matrix a41 (differentialProblem2.linearOperator ());
    dolfin::Vector a51 (differentialProblem2.rhs ());

    
    // ............................................................................................................ //
    dolfin::cout << dolfin::endl;
    dolfin::cout << "....................... THIRD CONSTRUCTOR ......................." << dolfin::endl;
    // define linear differential problem
    DCP::LinearDifferentialProblem<Poisson::BilinearForm, Poisson::LinearForm> differentialProblem3 (std::move(*mesh), std::move (*V));

//    dolfin::info (differentialProblem3.parameters, true);
    differentialProblem3.parameters ["desired_solver_type"] = "krylov_solver";
    differentialProblem3.parameters ["desired_solver_method"] = "gmres";
    differentialProblem3.parameters ["desired_solver_preconditioner"] = "ilu";
//    dolfin::info (differentialProblem3.parameters, true);
    
    differentialProblem3.update ();
//    dolfin::info (differentialProblem3.parameters, true);

    // set dirichlet bc
    differentialProblem3.addDirichletBC (dirichletBC);
    
    differentialProblem3.setCoefficient ("bilinear_form", c2, "c");
    differentialProblem3.setCoefficient ("linear_form", f2, "f");
    differentialProblem3.setCoefficient ("linear_form", g2, "g");
    
    differentialProblem3.setIntegrationSubdomains ("bilinear_form", 
                                                   dolfin::reference_to_no_delete_pointer (meshFacets),
                                                   DCP::SubdomainType::BOUNDARY_FACETS);
    differentialProblem3.setIntegrationSubdomains ("linear_form", 
                                                   dolfin::reference_to_no_delete_pointer (meshFacets),
                                                   DCP::SubdomainType::BOUNDARY_FACETS);
    
    // solve
    differentialProblem3.solve ();
    differentialProblem3.solve (true);
    differentialProblem3.solve ();
    differentialProblem3.parameters ["force_reassemble_system"] = true;
    differentialProblem3.solve ();
    differentialProblem3.parameters ["force_reassemble_system"] = false;
    differentialProblem3.solve ();
    differentialProblem3.parameters ["desired_solver_type"] = "lu_solver";
    differentialProblem3.parameters ["desired_solver_method"] = "default";
    differentialProblem3.solve ();
    dolfin::plot (differentialProblem3.solution ());
    
    dolfin::Form a12 (differentialProblem3.bilinearForm ());
    dolfin::Form a22 (differentialProblem3.linearForm ());
    dolfin::Matrix a42 (differentialProblem3.linearOperator ());
    dolfin::Vector a52 (differentialProblem3.rhs ());
    
    std::unique_ptr <DCP::AbstractDifferentialProblem> ptr (differentialProblem3.clone ());
    ptr->solve ();
    dolfin::plot (ptr->solution ());
     
    // =============================================================================================================== //
    // =============================================================================================================== //
    // =============================================================================================================== //
    dolfin::cout << dolfin::endl;   
    dolfin::cout << dolfin::endl;   
    dolfin::cout << dolfin::endl;   
    dolfin::cout << dolfin::endl;   
    dolfin::cout << dolfin::endl;   
    dolfin::cout << "************************************************************************************************" << dolfin::endl;
    dolfin::cout << "********************************* NONLINEAR PROBLEM TESTS **************************************" << dolfin::endl;
    dolfin::cout << "************************************************************************************************" << dolfin::endl;
    dolfin::cout << "------------------------------------------------------------------------------------------------" << dolfin::endl;
    dolfin::cout << "----------------------------------- NORMAL SOLUTION PROCESS ------------------------------------" << dolfin::endl;
    dolfin::cout << "------------------------------------------------------------------------------------------------" << dolfin::endl;
    // define mesh
    boost::shared_ptr<dolfin::UnitSquareMesh> NLmesh (new dolfin::UnitSquareMesh (20, 20));
    
    // define vector space
    boost::shared_ptr<NavierStokes::FunctionSpace> NLV (new NavierStokes::FunctionSpace (*NLmesh));
    
    // define variational form
    NavierStokes::ResidualForm NLF (*NLV);

    // define coefficients
    dolfin::Constant NLnu (1e-6);
    NLF.nu = NLnu;
    
    // define NLsolution
    dolfin::Function NLsolution (*NLV);
    
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
    boost::shared_ptr<dolfin::Function> NLinitialGuess (new dolfin::Function (NLsolution));
    NLF.set_coefficient ("trial", NLinitialGuess);
    
    NavierStokes::JacobianForm NLJ (*NLV, *NLV);
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
    
    dolfin::cout << dolfin::endl;   
    dolfin::cout << "------------------------------------------------------------------------------------------------" << dolfin::endl;
    dolfin::cout << "---------------------------- WITH NONLINEARDIFFERENTIALPROBLEM CLASS ---------------------------" << dolfin::endl;
    dolfin::cout << "------------------------------------------------------------------------------------------------" << dolfin::endl;
    dolfin::cout << dolfin::endl;
    dolfin::cout << "....................... FIRST CONSTRUCTOR ......................." << dolfin::endl;
    
    // define linear differential problem
    DCP::NonlinearDifferentialProblem<NavierStokes::ResidualForm, NavierStokes::JacobianForm> NLdifferentialProblem (*NLmesh, *NLV, "trial");

//  DCP::LinearDifferentialProblem<Poisson::BilinearForm, Poisson::LinearForm> differentialProblem
//      (std::move (mesh), std::move (V));

//    dolfin::info (NLdifferentialProblem.parameters, true);
    
    boost::shared_ptr<dolfin::Constant> NLnu2 (new dolfin::Constant (1e-6));
    NLdifferentialProblem.setCoefficient ("residual_form", NLnu2, "nu");
    NLdifferentialProblem.setCoefficient ("jacobian_form", NLnu2, "nu");
    
    for (auto i : NLboundaryConditions)
    {
        NLdifferentialProblem.addDirichletBC (*i);
    }
    
    NLdifferentialProblem.solve ();
    dolfin::plot (NLdifferentialProblem.solution ()[0]);
    dolfin::plot (NLdifferentialProblem.solution ()[0][0]);
    dolfin::plot (NLdifferentialProblem.solution ()[0][1]);
    dolfin::plot (NLdifferentialProblem.solution ()[1]);
    
    NLdifferentialProblem.parameters ("nonlinear_variational_solver")("newton_solver")["linear_solver"] = "cg";
//    dolfin::info (NLdifferentialProblem.parameters, true);
    
//    dolfin::list_linear_solver_methods ();
    NLdifferentialProblem.solve ();
    dolfin::plot (NLdifferentialProblem.solution ()[0]);
    dolfin::plot (NLdifferentialProblem.solution ()[0][0]);
    dolfin::plot (NLdifferentialProblem.solution ()[0][1]);
    dolfin::plot (NLdifferentialProblem.solution ()[1]);
    

    NLdifferentialProblem.parameters ("nonlinear_variational_solver")("newton_solver")["linear_solver"] = "default";
//    dolfin::info (NLdifferentialProblem.parameters, true);
    boost::shared_ptr<dolfin::Function> initGuess (new dolfin::Function (*NLV));
    NLdifferentialProblem.setCoefficient ("initial_guess", initGuess);
    NLdifferentialProblem.solve ();
    dolfin::plot (NLdifferentialProblem.solution ()[0]);
    dolfin::plot (NLdifferentialProblem.solution ()[0][0]);
    dolfin::plot (NLdifferentialProblem.solution ()[0][1]);
    dolfin::plot (NLdifferentialProblem.solution ()[1]);
    
    dolfin::cout << dolfin::endl;
    dolfin::cout << "....................... SECOND CONSTRUCTOR ......................." << dolfin::endl;
    
    // define linear differential problem
    DCP::NonlinearDifferentialProblem<NavierStokes::ResidualForm, NavierStokes::JacobianForm> NLdifferentialProblem2 (NLmesh, NLV, "trial");

//    dolfin::info (NLdifferentialProblem2.parameters, true);
    
    NLdifferentialProblem2.setCoefficient ("residual_form", NLnu2, "nu");
    NLdifferentialProblem2.setCoefficient ("jacobian_form", NLnu2, "nu");
    
    for (auto i : NLboundaryConditions)
    {
        NLdifferentialProblem2.addDirichletBC (*i);
    }
    
    NLdifferentialProblem2.solve ();
    dolfin::plot (NLdifferentialProblem2.solution ()[0]);
    dolfin::plot (NLdifferentialProblem2.solution ()[0][0]);
    dolfin::plot (NLdifferentialProblem2.solution ()[0][1]);
    dolfin::plot (NLdifferentialProblem2.solution ()[1]);
    
    NLdifferentialProblem2.parameters ("nonlinear_variational_solver")("newton_solver")["linear_solver"] = "cg";
//    dolfin::info (NLdifferentialProblem2.parameters, true);
    
//    dolfin::list_linear_solver_methods ();
    NLdifferentialProblem2.solve ();
    dolfin::plot (NLdifferentialProblem2.solution ()[0]);
    dolfin::plot (NLdifferentialProblem2.solution ()[0][0]);
    dolfin::plot (NLdifferentialProblem2.solution ()[0][1]);
    dolfin::plot (NLdifferentialProblem2.solution ()[1]);
    

    NLdifferentialProblem2.parameters ("nonlinear_variational_solver")("newton_solver")["linear_solver"] = "default";
//    dolfin::info (NLdifferentialProblem2.parameters, true);
    NLdifferentialProblem2.setCoefficient ("initial_guess", initGuess);
    NLdifferentialProblem2.solve ();
    dolfin::plot (NLdifferentialProblem2.solution ()[0]);
    dolfin::plot (NLdifferentialProblem2.solution ()[0][0]);
    dolfin::plot (NLdifferentialProblem2.solution ()[0][1]);
    dolfin::plot (NLdifferentialProblem2.solution ()[1]);
    
    dolfin::cout << dolfin::endl;
    dolfin::cout << "....................... THIRD CONSTRUCTOR ......................." << dolfin::endl;
    
    // define linear differential problem
    DCP::NonlinearDifferentialProblem<NavierStokes::ResidualForm, NavierStokes::JacobianForm> NLdifferentialProblem3 (std::move (*NLmesh), std::move (*NLV), "trial");

//    dolfin::info (NLdifferentialProblem3.parameters, true);
    
    NLdifferentialProblem3.setCoefficient ("residual_form", NLnu2, "nu");
    NLdifferentialProblem3.setCoefficient ("jacobian_form", NLnu2, "nu");
    
    for (auto i : NLboundaryConditions)
    {
        NLdifferentialProblem3.addDirichletBC (*i);
    }
    
    NLdifferentialProblem3.solve ();
    dolfin::plot (NLdifferentialProblem3.solution ()[0]);
    dolfin::plot (NLdifferentialProblem3.solution ()[0][0]);
    dolfin::plot (NLdifferentialProblem3.solution ()[0][1]);
    dolfin::plot (NLdifferentialProblem3.solution ()[1]);
    
    NLdifferentialProblem3.parameters ("nonlinear_variational_solver")("newton_solver")["linear_solver"] = "cg";
//    dolfin::info (NLdifferentialProblem3.parameters, true);
    
//    dolfin::list_linear_solver_methods ();
    NLdifferentialProblem3.solve ();
    dolfin::plot (NLdifferentialProblem3.solution ()[0]);
    dolfin::plot (NLdifferentialProblem3.solution ()[0][0]);
    dolfin::plot (NLdifferentialProblem3.solution ()[0][1]);
    dolfin::plot (NLdifferentialProblem3.solution ()[1]);
    

    NLdifferentialProblem3.parameters ("nonlinear_variational_solver")("newton_solver")["linear_solver"] = "default";
//    dolfin::info (NLdifferentialProblem3.parameters, true);
    NLdifferentialProblem3.setCoefficient ("initial_guess", initGuess);
    NLdifferentialProblem3.solve ();
    dolfin::plot (NLdifferentialProblem3.solution ()[0]);
    dolfin::plot (NLdifferentialProblem3.solution ()[0][0]);
    dolfin::plot (NLdifferentialProblem3.solution ()[0][1]);
    dolfin::plot (NLdifferentialProblem3.solution ()[1]);
    
    std::unique_ptr <DCP::AbstractDifferentialProblem> NLptr (NLdifferentialProblem3.clone ());
    NLptr -> solve ();
    dolfin::plot (NLptr->solution ()[0]);
    
    DCP::NonlinearDifferentialProblem<NavierStokes::ResidualForm, NavierStokes::JacobianForm> NLdifferentialProblem4 = NLdifferentialProblem2;

    std::unique_ptr <DCP::AbstractDifferentialProblem> Lptr1 (differentialProblem.clone ());
    
    NLdifferentialProblem3.solve ();
    dolfin::plot (NLdifferentialProblem3.solution ()[0]);
    dolfin::plot (NLdifferentialProblem3.solution ()[0][0]);
    dolfin::plot (NLdifferentialProblem3.solution ()[0][1]);
    dolfin::plot (NLdifferentialProblem3.solution ()[1]);
    
    differentialProblem.solve ();
    dolfin::plot (differentialProblem.solution ());
    
    std::map <std::string, dolfin::DirichletBC> mappa = differentialProblem.dirichletBCs ();
    for (auto &i : mappa)
        std::cout << i.first << std::endl;
    
    mappa = NLdifferentialProblem3.dirichletBCs ();
    for (auto &i : mappa)
        std::cout << i.first << std::endl;
    
    differentialProblem.removeDirichletBC ("culo");
    auto res = differentialProblem.removeDirichletBC ("my_name");
    std::cout << "res = " << res << std::endl;
    
    
    return 0;
}
