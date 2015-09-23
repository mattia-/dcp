// Begin demo

#include <dolfin.h>
#include "myNavierstokesTimeCurv.h"
#include "computeFreeSurfaceStress_onlyTP.h"
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
  dolfin::RectangleMesh mesh(0,0, lx,ly, nx, ny);

  // Time stepping and model parameters
//  dolfin::Constant dt(ustar*2.5*10.0/nx);//0.05);
#if defined(GerbeauLelievre) && !defined(SprittlesShikhmurzaev) && !defined(Yamamoto)
  dolfin::Constant dt(0.1);//0.05);
  dolfin::Constant nu (1.95 / 0.81);
std::cerr << "Reynolds number = " << re << std::endl;
  dolfin::Constant gamma(5.5);
//  dolfin::Constant w(0,0);
  dolfin::Constant beta(36);
  dolfin::Constant cosThetaS (cos(3.14159265/2.0));
  dolfin::Constant stressBelow(0,0);
#elif defined(Yamamoto) && !defined(GerbeauLelievre) && !defined(SprittlesShikhmurzaev)
  dolfin::Constant dt(0.00002);//0.05);
  #ifdef vanMourik
      dolfin::Constant nu (2.91e-6);
    std::cerr << "Reynolds number = " << ustar*lx/2.91e-6 << std::endl;
      dolfin::Constant gamma(2.06e-5);
      dolfin::Constant cosThetaS (cos(53.6*3.14159265/180.0));
      dolfin::Constant stressBelow(0,0);
  #else
  dolfin::Constant nu (2.081e-2 / rho);
std::cerr << "Reynolds number = " << ustar*lx*2.081e-2/rho << std::endl;
  dolfin::Constant gamma(4.36e-5);
  dolfin::Constant cosThetaS (cos(69.8*3.14159265/180.0));
  dolfin::Constant stressBelow(0,rho*9.81*ly);
  #endif
//  dolfin::Constant w(0,0);
  dolfin::Constant beta(36);
#endif
  dolfin::Constant t_partialOmega(0,1);
  #ifdef vanMourik
      dolfin::Constant gravityVersor(0,0);
  #else
  dolfin::Constant gravityVersor(0,-1);
  #endif
  dcp::TimeDependentExpression stressAbove(2,(BoundaryStressEvaluator (0,0)));

  double t0 = 0.0;
  double T  = 50;

  // Create user-defined nonlinear problem
  myNavierstokesTimeCurv::FunctionSpace V (mesh);
  Ivan::UflToNewton < myNavierstokesTimeCurv::FunctionSpace, computeFreeSurfaceStress_onlyTP::FunctionSpace,
                      myNavierstokesTimeCurv::ResidualForm, myNavierstokesTimeCurv::JacobianForm,
                      computeFreeSurfaceStress_onlyTP::LinearForm >
       timeSteppingProblem(V,"trial");//,"dtrial");

  // Create nonlinear solver, set parameters and assign it to timeSteppingProblem
  dolfin::NewtonSolver newton_solver;
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
 // LateralInterior wallBoundaryInterior;
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
  navierStokesProblem.addDirichletBC (
  #ifdef vanMourik
            dolfin::DirichletBC (* (* navierStokesProblem.functionSpace())[0], noSlipDirichletBC, inflowBoundary),
  #else
        dolfin::DirichletBC (* (* (* navierStokesProblem.functionSpace())[0])[0], freeSlipDirichletBC, inflowBoundary),
  #endif
        std::string("cut"));

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
  navierStokesProblem.setCoefficient ("residual_form", dolfin::reference_to_no_delete_pointer (nu),"nu");
//  navierStokesProblem.setCoefficient ("residual_form", dolfin::reference_to_no_delete_pointer (dt),"dt");
  navierStokesProblem.setCoefficient ("residual_form", dolfin::reference_to_no_delete_pointer (gamma),"gamma");
//  navierStokesProblem.setCoefficient ("residual_form", dolfin::reference_to_no_delete_pointer (nonlinearProblem.previousSolution()[0]),"u_old");
//  navierStokesProblem.setCoefficient ("residual_form", dolfin::reference_to_no_delete_pointer (nonlinearProblem.solution()),"trial");
  navierStokesProblem.setCoefficient ("residual_form", dolfin::reference_to_no_delete_pointer (beta),"beta");
  navierStokesProblem.setCoefficient ("residual_form", dolfin::reference_to_no_delete_pointer (stressAbove),"stressAbove");
  navierStokesProblem.setCoefficient ("residual_form", dolfin::reference_to_no_delete_pointer (stressBelow),"stressBelow");
  navierStokesProblem.setCoefficient ("residual_form", dolfin::reference_to_no_delete_pointer (*(new WallVelocity)),"wallVelocity");
  navierStokesProblem.setCoefficient ("residual_form", dolfin::reference_to_no_delete_pointer (gravityVersor),"gravityVersor");
  navierStokesProblem.setCoefficient ("jacobian_form", dolfin::reference_to_no_delete_pointer (nu),"nu");
//  navierStokesProblem.setCoefficient ("jacobian_form", dolfin::reference_to_no_delete_pointer (dt),"dt");
//  navierStokesProblem.setCoefficient ("jacobian_form", dolfin::reference_to_no_delete_pointer (nonlinearProblem.solution()),"trial");
//  navierStokesProblem.setCoefficient ("additional_form", dolfin::reference_to_no_delete_pointer (dt),"dt");
  navierStokesProblem.setCoefficient ("jacobian_form", dolfin::reference_to_no_delete_pointer (beta),"beta");
  navierStokesProblem.setCoefficient ("additional_form", dolfin::reference_to_no_delete_pointer (gamma),"gamma");
  navierStokesProblem.setCoefficient ("additional_form", dolfin::reference_to_no_delete_pointer (cosThetaS),"cosThetaS");
  navierStokesProblem.setCoefficient ("additional_form", dolfin::reference_to_no_delete_pointer (t_partialOmega),"t_partialOmega");

  // Solution functions
//  dolfin::Function& previousSolution = nonlinearProblem.solution();

  // Save initial condition to file
  dolfin::Constant initialCondition (0,0,0);
  navierStokesProblem.setInitialSolution (initialCondition);
  dolfin::File file("solution.pvd");
  const dolfin::Function& solution = navierStokesProblem.solution();
  file << std::make_pair(&solution, t0);

  // Solve
  navierStokesProblem.solve ();

  // Plot and save solution
  const dolfin::Function& final_solution = navierStokesProblem.solution();
  dolfin::plot(final_solution[0]);
  dolfin::interactive();
  file << std::make_pair(&final_solution, T);

  // Post-processing [GerbeauLelievre]


  return 0;
}

