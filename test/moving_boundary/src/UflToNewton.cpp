// Begin demo

#include <dolfin.h>
#include "myNavierstokesTimeCurv.h"
#include "computeFreeSurfaceStress_onlyTP.h"
#include "UflToNewton.h"
#include "utilities.h"

/*// Initial conditions
class InitialConditions : public dolfin::Expression
{
public:

  InitialConditions() : dolfin::Expression(2)
  {
    dolfin::seed(2 + dolfin::MPI::rank(MPI_COMM_WORLD));
  }

  void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x) const
  {
    values[0]= 1.0;//0.63 + 0.02*(0.5 - dolfin::rand());
    values[1]= 0.0;
  }

};*/

// Subdomains and expressions for BCs

int main(int argc, char* argv[])
{
  dolfin::init(argc, argv);

  // Mesh
//  dolfin::UnitSquareMesh mesh(20, 20);
  dolfin::RectangleMesh mesh(0,0, 1,1, 2, 2);

  // Time stepping and model parameters
  dolfin::Constant dt(0.1);
  dolfin::Constant nu(1.0e-01);
  dolfin::Constant gamma(7.3e-5);
  dolfin::Constant w(0,0);

  double t = 0.0;
  double T = 10*dt;

  // Create user-defined nonlinear problem
  Ivan::UflToNewton < myNavierstokesTimeCurv::FunctionSpace, computeFreeSurfaceStress_onlyTP::FunctionSpace,
                      myNavierstokesTimeCurv::ResidualForm, myNavierstokesTimeCurv::JacobianForm,
                      computeFreeSurfaceStress_onlyTP::LinearForm >
       nonlinearProblem(mesh);
std::cerr << "OK fino a " << __LINE__ << std::endl;
  nonlinearProblem.setCoefficient ("residual_form", dolfin::reference_to_no_delete_pointer (nu),"nu");
  nonlinearProblem.setCoefficient ("residual_form", dolfin::reference_to_no_delete_pointer (dt),"dt");
  nonlinearProblem.setCoefficient ("residual_form", dolfin::reference_to_no_delete_pointer (gamma),"gamma");
  nonlinearProblem.setCoefficient ("residual_form", dolfin::reference_to_no_delete_pointer (nonlinearProblem.previousSolution()[0]),"u_old");
  nonlinearProblem.setCoefficient ("residual_form", dolfin::reference_to_no_delete_pointer (nonlinearProblem.solution()),"trial");
  nonlinearProblem.setCoefficient ("jacobian_form", dolfin::reference_to_no_delete_pointer (nu),"nu");
  nonlinearProblem.setCoefficient ("jacobian_form", dolfin::reference_to_no_delete_pointer (dt),"dt");
  nonlinearProblem.setCoefficient ("jacobian_form", dolfin::reference_to_no_delete_pointer (nonlinearProblem.solution()),"trial");
  nonlinearProblem.setCoefficient ("additional_form", dolfin::reference_to_no_delete_pointer (dt),"dt");
  nonlinearProblem.setCoefficient ("additional_form", dolfin::reference_to_no_delete_pointer (gamma),"gamma");

  // Boundary conditions
  dolfin::Constant noSlipDirichletBC (0.0, 0.0);
  dolfin::Constant freeSlipDirichletBC (0.0);
  InflowDirichletData inflowDirichletBC;

  BottomBoundary inflowBoundary;
  LateralInterior wallBoundary;
  TopBoundary freeSurface;
  TriplePointLeft triplePointLeft;
  TriplePointRight triplePointRight;
  TriplePointLeftVertex triplePointLeftVertex;
  TriplePointRightVertex triplePointRightVertex;

  nonlinearProblem.addDirichletBC (std::string("inflow"),
        dolfin::DirichletBC (* nonlinearProblem.functionSpace()[0], inflowDirichletBC, inflowBoundary));
  nonlinearProblem.addDirichletBC (std::string("wall"),
        dolfin::DirichletBC (* (* nonlinearProblem.functionSpace()[0])[0], freeSlipDirichletBC, wallBoundary));
//  nonlinearProblem.addDirichletBC (std::string("moving"), dolfin::DirichletBC (* nonlinearProblem.functionSpace()[0], movingLidVelocity, movingLid));
//  nonlinearProblem.addDirichletBC (std::string("fixed"), dolfin::DirichletBC (* nonlinearProblem.functionSpace()[0], noSlipDirichletBC, fixedWalls));

  dolfin::FacetFunction<std::size_t> meshFacets (mesh);
  meshFacets.set_all (0);
  freeSurface.mark (meshFacets, 1);
  dolfin::FacetFunction<std::size_t> additionalMeshFacets (mesh);
  additionalMeshFacets.set_all (0);
  triplePointLeft.mark(additionalMeshFacets,4);
  triplePointRight.mark(additionalMeshFacets,5);
  dolfin::VertexFunction<std::size_t> additionalMeshVertices (mesh);
  additionalMeshVertices.set_all (0);
  triplePointLeftVertex.mark(additionalMeshVertices,40);
  triplePointRightVertex.mark(additionalMeshVertices,50);
std::cerr << "OK fino a " << __LINE__ << std::endl;
  nonlinearProblem.setIntegrationSubdomain ("residual_form",
        dolfin::reference_to_no_delete_pointer (meshFacets), dcp::SubdomainType::BOUNDARY_FACETS);
  nonlinearProblem.setIntegrationSubdomain ("additional_form",
        dolfin::reference_to_no_delete_pointer (additionalMeshFacets), dcp::SubdomainType::BOUNDARY_FACETS);
  nonlinearProblem.setIntegrationSubdomain ("additional_form",
        dolfin::reference_to_no_delete_pointer (additionalMeshVertices), dcp::SubdomainType::VERTICES);
std::cerr << "OK fino a " << __LINE__ << std::endl;
dolfin::plot(meshFacets, "meshFacets");dolfin::plot(additionalMeshFacets, "additionalMeshFacets");//dolfin::interactive();

  // Solution functions
  dolfin::Function& solution = nonlinearProblem.solution();
  dolfin::Function& previousSolution = nonlinearProblem.solution();

  // Create nonlinear solver and set parameters
  dolfin::NewtonSolver newton_solver;
  newton_solver.parameters["linear_solver"] = "lu";
  newton_solver.parameters["convergence_criterion"] = "incremental";
  newton_solver.parameters["maximum_iterations"] = 50;
/*  newton_solver.parameters["relative_tolerance"] = 1e-6;
  newton_solver.parameters["absolute_tolerance"] = 1e-15;*/

  // Save initial condition to file
  dolfin::Constant initialCondition (0,0,0);
  previousSolution = ( solution = initialCondition );
  dolfin::File file("cahn_hilliard.pvd");
  file << std::make_pair(&solution, t);

  // Solve
  while (t < T)
  {
    // Update for next time step
    t += dt;
    *previousSolution.vector() = *solution.vector();

    // Updating coefficients
    nonlinearProblem.setCoefficient ("residual_form", dolfin::reference_to_no_delete_pointer (w), "w");
    nonlinearProblem.setCoefficient ("residual_form", dolfin::reference_to_no_delete_pointer (previousSolution[0]),"u_old");
    nonlinearProblem.setCoefficient ("jacobian_form", dolfin::reference_to_no_delete_pointer (w), "w");

    // Solve
    newton_solver.solve(nonlinearProblem, *solution.vector());

    // Save function to file
    file << std::pair<const dolfin::Function*, double>(&solution, t);
  }

  // Plot solution
  dolfin::plot(solution[0]);
//  dolfin::plot(solution[1]);
  dolfin::interactive();

  return 0;
}

