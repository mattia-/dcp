// Begin demo

#include <dolfin.h>
#include "myNavierstokesTimeCurv_dimensionless.h"
#include "computeFreeSurfaceStress_onlyTP_dimensionless.h"
#include "UflToNewton.h"
#include "utilities.h"
#include "MovingTimeDependentProblem.h"
#include <dcp/subdomains/Subdomain.h>
#include <math.h>

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
  dolfin::RectangleMesh mesh(0,0, lx,ly, nx, yxratio*nx,"crossed");
dolfin::ALE meshMover;
InitialDisplacement initialDisplacementExpr;
//meshMover.move(mesh,initialDisplacementExpr);

  // Time stepping and model parameters
  dolfin::Constant dt(0.001);//0.1*ustar*2.5*10.0/nx);//0.05);
std::cerr << "Reynolds number = " << re << std::endl;
  dolfin::Constant reynoldsNr (re);
  dolfin::Constant stokesNr (st);
  dolfin::Constant capillaryNr (ca);
  dolfin::Constant beta(1.0e4/(re*re));
  dolfin::Constant cosThetaS (cos(3.14159265/3.0));
  dolfin::Constant t_partialOmega(0,1);
  dolfin::Constant gravityVersor(0,0);//-1);
  dolfin::Constant stressBelow(0,0);//ly);
  dcp::TimeDependentExpression stressAbove(2,(BoundaryStressEvaluator (0,0)));//1)));
  dolfin::Constant wallVelocity (0,-1);

  double t0 = 0.0;
  double T  = t0+750*dt;

  // Create user-defined nonlinear problem
  myNavierstokesTimeCurv_dimensionless::FunctionSpace V (mesh);
  Ivan::UflToNewton < myNavierstokesTimeCurv_dimensionless::FunctionSpace, computeFreeSurfaceStress_onlyTP_dimensionless::FunctionSpace,
                      myNavierstokesTimeCurv_dimensionless::ResidualForm, myNavierstokesTimeCurv_dimensionless::JacobianForm,
                      computeFreeSurfaceStress_onlyTP_dimensionless::LinearForm >
       timeSteppingProblem(V,"trial");//,"dtrial");

  // Create nonlinear solver, set parameters and assign it to timeSteppingProblem
  dolfin::NewtonSolver newton_solver;
  newton_solver.parameters["linear_solver"] = "lu";
  newton_solver.parameters["convergence_criterion"] = "incremental";
  newton_solver.parameters["maximum_iterations"] = 50;
/*  newton_solver.parameters["relative_tolerance"] = 1e-6;
  newton_solver.parameters["absolute_tolerance"] = 1e-15;*/
  timeSteppingProblem.setNonlinearSolver (dolfin::reference_to_no_delete_pointer(newton_solver));


  Ivan::MovingTimeDependentProblem navierStokesProblem (dolfin::reference_to_no_delete_pointer (timeSteppingProblem),
                                                   t0,
                                                   dt, 
                                                   T, 
                                                   {"residual_form", "jacobian_form", "additional_form"},
                                                   {"residual_form"}
//current//                                                  {"residual_form", "jacobian_form"}
                                                  );
         
std::cerr << "OK fino a " << __LINE__ << std::endl;

  // Boundary conditions
  dolfin::Constant noSlipDirichletBC (0.0, 0.0);
  dolfin::Constant freeSlipDirichletBC (0.0);
//inflowTime//  InflowDirichletData inflowDirichletBC;

//inflowTime//  BottomBoundary inflowBoundary;
  LateralWall wallBoundary;
  LateralInterior wallBoundaryInterior;
  TopBoundary freeSurface;
  TriplePointLeft triplePointLeft;
  TriplePointRight triplePointRight;
  TriplePointLeftVertex triplePointLeftVertex;
  TriplePointRightVertex triplePointRightVertex;

//inflowTime//  navierStokesProblem.addDirichletBC (
//inflowTime//        dolfin::DirichletBC (* (* navierStokesProblem.functionSpace())[0], inflowDirichletBC, inflowBoundary),
//inflowTime//        std::string("inflow"));
  dcp::Subdomain inflowBoundary ((BottomBoundaryEvaluator ()));
//noInflow//  dcp::TimeDependentExpression inflowDirichletBC(2,(InflowDirichletBCEvaluator ()));
//noInflow//  navierStokesProblem.addTimeDependentDirichletBC (inflowDirichletBC, inflowBoundary, 0, "inflowBC");

  navierStokesProblem.addDirichletBC (
        dolfin::DirichletBC (* (* (* navierStokesProblem.functionSpace())[0])[0], freeSlipDirichletBC, wallBoundary),
//        dolfin::DirichletBC (* (* navierStokesProblem.functionSpace())[0], noSlipDirichletBC, wallBoundary),
        std::string("wall"));

  dolfin::FacetFunction<std::size_t> meshFacets (navierStokesProblem.mesh());
  meshFacets.set_all (0);
  freeSurface.mark (meshFacets, 1);
  wallBoundary.mark (meshFacets, 2);
  inflowBoundary.mark (meshFacets, 3);
//!!! NB uso mesh invece di navierStokesProblem.mesh() perche' altrimenti l'assembler della additionalForm_ da' un errore di inconsistenza tra
//       la mesh di additionalFunctionSpace_ e quella delle misure (additionalMeshFacets, additionalMeshVertices).
//       Perche'? boh...
  dolfin::FacetFunction<std::size_t> additionalMeshFacets (mesh);
  additionalMeshFacets.set_all (0);
  triplePointLeft.mark(additionalMeshFacets,4);
  triplePointRight.mark(additionalMeshFacets,5);
/*addVert//  dolfin::VertexFunction<std::size_t> additionalMeshVertices (mesh);
  additionalMeshVertices.set_all (0);
  triplePointLeftVertex.mark(additionalMeshVertices,40);
  triplePointRightVertex.mark(additionalMeshVertices,50);
*///addVert
std::cerr << "OK fino a " << __LINE__ << std::endl;
  navierStokesProblem.setIntegrationSubdomain ("residual_form",
        dolfin::reference_to_no_delete_pointer (meshFacets), dcp::SubdomainType::BOUNDARY_FACETS);
  navierStokesProblem.setIntegrationSubdomain ("additional_form",
        dolfin::reference_to_no_delete_pointer (additionalMeshFacets), dcp::SubdomainType::BOUNDARY_FACETS);
/*addVert//  navierStokesProblem.setIntegrationSubdomain ("additional_form",
        dolfin::reference_to_no_delete_pointer (additionalMeshVertices), dcp::SubdomainType::VERTICES);
*///addVert
std::cerr << "OK fino a " << __LINE__ << std::endl;
dolfin::plot(meshFacets, "meshFacets");dolfin::plot(additionalMeshFacets, "additionalMeshFacets");dolfin::interactive();

  // Coefficients
  navierStokesProblem.setCoefficient ("residual_form", dolfin::reference_to_no_delete_pointer (reynoldsNr),"reynoldsNr");
  navierStokesProblem.setCoefficient ("residual_form", dolfin::reference_to_no_delete_pointer (stokesNr),"stokesNr");
  navierStokesProblem.setCoefficient ("residual_form", dolfin::reference_to_no_delete_pointer (capillaryNr),"capillaryNr");
  navierStokesProblem.setCoefficient ("residual_form", dolfin::reference_to_no_delete_pointer (beta),"beta");
  navierStokesProblem.setCoefficient ("residual_form", dolfin::reference_to_no_delete_pointer (stressAbove),"stressAbove");
  navierStokesProblem.setCoefficient ("residual_form", dolfin::reference_to_no_delete_pointer (stressBelow),"stressBelow");
  navierStokesProblem.setCoefficient ("residual_form", dolfin::reference_to_no_delete_pointer (gravityVersor),"gravityVersor");
  navierStokesProblem.setCoefficient ("residual_form", dolfin::reference_to_no_delete_pointer (wallVelocity),"wallVelocity");
  navierStokesProblem.setCoefficient ("jacobian_form", dolfin::reference_to_no_delete_pointer (reynoldsNr),"reynoldsNr");
  navierStokesProblem.setCoefficient ("jacobian_form", dolfin::reference_to_no_delete_pointer (beta),"beta");
  navierStokesProblem.setCoefficient ("additional_form", dolfin::reference_to_no_delete_pointer (reynoldsNr),"reynoldsNr");
  navierStokesProblem.setCoefficient ("additional_form", dolfin::reference_to_no_delete_pointer (capillaryNr),"capillaryNr");
  navierStokesProblem.setCoefficient ("additional_form", dolfin::reference_to_no_delete_pointer (cosThetaS),"cosThetaS");
  navierStokesProblem.setCoefficient ("additional_form", dolfin::reference_to_no_delete_pointer (t_partialOmega),"t_partialOmega");

  // Solution functions
//  dolfin::Function& previousSolution = nonlinearProblem.solution();

  // Save initial condition to file
  dolfin::Constant initialCondition (0,0,0);
  dolfin::Function initialSolution (navierStokesProblem.functionSpace());
  initialSolution = initialCondition;
  navierStokesProblem.setInitialSolution (initialSolution);
  dolfin::File file("solution.pvd");
  const dolfin::Function& solution = navierStokesProblem.solution();
  file << std::make_pair(&solution, t0);

  // Solve
  navierStokesProblem.step ();
  navierStokesProblem.clear ();
  navierStokesProblem.solve ();

  // Plot and save solution
  const dolfin::Function& final_solution = navierStokesProblem.solution();
  dolfin::plot(final_solution[0]);
  dolfin::interactive();
  file << std::make_pair(&final_solution, T);

  return 0;
}

