#include <iostream>
#include <dolfin.h>
#include "primal.h"
#include "adjoint.h"
#include <DifferentialProblem/NonlinearDifferentialProblem.hpp>
#include <DifferentialProblem/LinearDifferentialProblem.hpp>
#include <DifferentialProblem/CompositeDifferentialProblem.hpp>

namespace primal
{
    class InflowBoundary : public dolfin::SubDomain
    {
        bool inside (const dolfin::Array<double>& x, bool on_boundary) const
        {
            return x[0] < (0 + DOLFIN_EPS) && on_boundary;
        }
    };

    class GammaSD : public dolfin::SubDomain
    {
        bool inside (const dolfin::Array<double>& x, bool on_boundary) const
        {
            return (x[1] < (0 + DOLFIN_EPS) || x[1] > (2 - DOLFIN_EPS)) && on_boundary;
        }
    };
}

namespace adjoint
{
    class DirichletBoundary : public dolfin::SubDomain
    {
        bool inside (const dolfin::Array<double>& x, bool on_boundary) const
        {
            return x[0] < (0 + DOLFIN_EPS) || x[1] < (0 + DOLFIN_EPS) || x[1] > (2 - DOLFIN_EPS);
        }
    };

    class RobinBoundary : public dolfin::SubDomain
    {
        bool inside (const dolfin::Array<double>& x, bool on_boundary) const
        {
            return x[0] > (5 - DOLFIN_EPS);
        }
    };
}

int main ()
{
    dolfin::set_log_level (dolfin::DBG);

    // mesh
    dolfin::RectangleMesh mesh (0.0, 0.0, 5.0, 2.0, 50, 20);

    dolfin::plot (mesh);

    // constant
    dolfin::Constant nu (1e-5);

    // define primal problem
    primal::FunctionSpace V (mesh);
    controlproblem::NonlinearDifferentialProblem <primal::ResidualForm, primal::JacobianForm> 
        primalProblem (dolfin::reference_to_no_delete_pointer (mesh), 
                       dolfin::reference_to_no_delete_pointer (V),
                       "trial");

    primalProblem.setCoefficient ("residual_form", dolfin::reference_to_no_delete_pointer (nu), "nu");
    primalProblem.setCoefficient ("jacobian_form", dolfin::reference_to_no_delete_pointer (nu), "nu");

    dolfin::Constant primal_inflowDirichletBC (1.0, 0.0);
    dolfin::Constant primal_symmetryDirichletBC (0.0);

    primal::InflowBoundary primal_inflowBoundary;
    primal::GammaSD primal_gammaSD;
    primalProblem.addDirichletBC (dolfin::DirichletBC (*V[0], primal_inflowDirichletBC, primal_inflowBoundary));
    primalProblem.addDirichletBC (dolfin::DirichletBC (*(*V[0])[1], primal_symmetryDirichletBC, primal_gammaSD));

    primalProblem.solve ();

    dolfin::plot (primalProblem.solution ()[0]);
    dolfin::plot (primalProblem.solution ()[1]);

    dolfin::UnitSquareMesh mesh2 (10, 10);
    primal::FunctionSpace V2 (mesh2);
    dolfin::Function u (V2[0]->collapse ());
    u.interpolate (primalProblem.solution () [0]);
    dolfin::plot (u);

    controlproblem::LinearDifferentialProblem <adjoint::BilinearForm, adjoint::LinearForm> 
        adjointProblem (dolfin::reference_to_no_delete_pointer (mesh),
                        dolfin::reference_to_no_delete_pointer (V));

    adjointProblem.setCoefficient ("bilinear_form", dolfin::reference_to_no_delete_pointer (nu), "nu");
    adjointProblem.setCoefficient ("bilinear_form", dolfin::reference_to_no_delete_pointer (primalProblem.solution ()[0]), "u");
    adjointProblem.setCoefficient ("linear_form", boost::shared_ptr<dolfin::Constant> (new dolfin::Constant (1.0, 0.0)), "F_u");
    adjointProblem.setCoefficient ("linear_form", boost::shared_ptr<dolfin::Constant> (new dolfin::Constant (1.0)), "F_p");

    dolfin::Constant adjoint_dirichletBC (0.0, 0.0);
    adjoint::DirichletBoundary adjoint_dirichletBoundary;
    adjointProblem.addDirichletBC (dolfin::DirichletBC (*V[0], adjoint_dirichletBC, adjoint_dirichletBoundary));

    adjoint::RobinBoundary adjoint_robinBoundary;
    dolfin::FacetFunction<std::size_t> adjoint_meshFacets (mesh);
    adjoint_meshFacets.set_all (1);
    adjoint_robinBoundary.mark (adjoint_meshFacets, 0);

    adjointProblem.setIntegrationSubdomains ("bilinear_form", 
                                             dolfin::reference_to_no_delete_pointer (adjoint_meshFacets),
                                             controlproblem::SubdomainType::BOUNDARY_FACETS);

    adjointProblem.solve ();

    dolfin::plot (adjointProblem.solution () [0]);
    dolfin::plot (adjointProblem.solution () [1]);

    // ============================================================================================================ //
    // ============================================================================================================ //
    // ============================================================================================================ //
    // ============================================================================================================ //
    // ============================================================================================================ //
    // ============================================================================================================ //

    // define problems
    controlproblem::NonlinearDifferentialProblem <primal::ResidualForm, primal::JacobianForm> 
        primalProblem2 (dolfin::reference_to_no_delete_pointer (mesh), 
                        dolfin::reference_to_no_delete_pointer (V),
                        "trial");

    controlproblem::LinearDifferentialProblem <adjoint::BilinearForm, adjoint::LinearForm> 
        adjointProblem2 (dolfin::reference_to_no_delete_pointer (mesh),
                         dolfin::reference_to_no_delete_pointer (V));

    // creating composite differential problem
    controlproblem::CompositeDifferentialProblem problems;

    problems.addProblem ("primal", primalProblem2);
    problems.addProblem ("adjoint", adjointProblem2);

    problems["primal"].setCoefficient ("residual_form", dolfin::reference_to_no_delete_pointer (nu), "nu");
    problems["primal"].setCoefficient ("jacobian_form", dolfin::reference_to_no_delete_pointer (nu), "nu");

    problems["primal"].addDirichletBC (dolfin::DirichletBC (*V[0], primal_inflowDirichletBC, primal_inflowBoundary));
    problems["primal"].addDirichletBC (dolfin::DirichletBC (*(*V[0])[1], primal_symmetryDirichletBC, primal_gammaSD));


    problems["adjoint"].setCoefficient ("bilinear_form", dolfin::reference_to_no_delete_pointer (nu), "nu");
    problems["adjoint"].setCoefficient ("linear_form", boost::shared_ptr<dolfin::Constant> (new dolfin::Constant (1.0, 0.0)), "F_u");
    problems["adjoint"].setCoefficient ("linear_form", boost::shared_ptr<dolfin::Constant> (new dolfin::Constant (1.0)), "F_p");

    problems.linkProblems ("adjoint", "u", "bilinear_form", "primal", 0);
    problems.print ();
    
    problems["adjoint"].addDirichletBC (dolfin::DirichletBC (*V[0], adjoint_dirichletBC, adjoint_dirichletBoundary));

    problems["adjoint"].setIntegrationSubdomains ("bilinear_form", 
                                                  dolfin::reference_to_no_delete_pointer (adjoint_meshFacets),
                                                  controlproblem::SubdomainType::BOUNDARY_FACETS);

    problems.solve ("primal");
    problems.solve ("adjoint");
    
    problems.solve ();
    problems.solve ();


    dolfin::plot (problems.solution ("primal")[0]);
    dolfin::plot (problems.solution ("primal")[1]);
    dolfin::plot (problems.solution ("adjoint")[0]);
    dolfin::plot (problems.solution ("adjoint")[1]);

    dolfin::interactive ();

    return 0;
}
