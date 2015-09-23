// Begin demo

#include <dolfin.h>
#include "myNavierstokesTimeCurvLinear.h"
#include "computeFreeSurfaceStress_onlyTP.h"
#include "UflToNewtonLinear.h"
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

template <typename T>
class MyFacetFunction : public dolfin::FacetFunction<T>
{
    using dolfin::FacetFunction<T>::FacetFunction;

  public:
    std::shared_ptr<const dolfin::Mesh> mesh()
    {
      dolfin_assert(this->_mesh);
      this->_mesh.reset (pb_->mesh());
      dolfin_assert(this->_mesh);
      return this->_mesh;
    }

    void setPb (Ivan::MovingTimeDependentProblem& pb)
    {
      pb_.reset (&pb);
    }

  private:
    std::shared_ptr<Ivan::MovingTimeDependentProblem> pb_;
};

// Subdomains and expressions for BCs

int main(int argc, char* argv[])
{
  dolfin::init(argc, argv);
dolfin::parameters["allow_extrapolation"] = true;

  // Mesh
//  dolfin::UnitSquareMesh mesh(20, 20);
  dolfin::RectangleMesh mesh(0,0, lx,ly, nx, ny);

InitialDisplacement initialDisplacement;
dolfin::plot(initialDisplacement,mesh); dolfin::interactive();
mesh.move(initialDisplacement);
dolfin::plot((XY()),mesh); dolfin::interactive();

  // Time stepping and model parameters
//  dolfin::Constant dt(ustar*2.5*10.0/nx);//0.05);
//  dolfin::Constant w(0,0);
  dolfin::Constant beta(betaVal);
  dolfin::Constant cosThetaS (cosThetaSVal);
#if defined(GerbeauLelievre) && !defined(SprittlesShikhmurzaev) && !defined(Yamamoto)
  dolfin::Constant dt(0.1);//0.05);
  dolfin::Constant nu (1.95 / 0.81);
std::cerr << "Reynolds number = " << re << std::endl;
  dolfin::Constant gamma(5.5);
  dolfin::Constant stressBelow(0,0);
#elif defined(Yamamoto) && !defined(GerbeauLelievre) && !defined(SprittlesShikhmurzaev)
  dolfin::Constant dt(0.00002);//0.05);
  #ifdef vanMourik
      dolfin::Constant nu (2.91e-6);
    std::cerr << "Reynolds number = " << ustar*lx/2.91e-6 << std::endl;
      dolfin::Constant gamma(2.06e-5);
      dolfin::Constant stressBelow(0,0);
  #else
  dolfin::Constant nu (2.081e-2 / rho);
std::cerr << "Reynolds number = " << ustar*lx*rho/2.081e-2 << std::endl;
  dolfin::Constant gamma(4.36e-5);
  dolfin::Constant stressBelow(0,9.81*ly);
  #endif
#endif
  dolfin::Constant t_partialOmega(0,1);
  #ifdef vanMourik
      dolfin::Constant gravityVersor(0,0);
  #else
  dolfin::Constant gravityVersor(0,-1);
  #endif
  dcp::TimeDependentExpression stressAbove(2,(BoundaryStressEvaluator (0,0)));

  double t0 = 0.0;
  double T  = t0+200*dt;
  //double T  = t0+10000*dt;

  // Create user-defined linearized problem
  myNavierstokesTimeCurvLinear::FunctionSpace V (mesh);
  Ivan::UflToNewtonLinear < myNavierstokesTimeCurvLinear::FunctionSpace, computeFreeSurfaceStress_onlyTP::FunctionSpace,
                      myNavierstokesTimeCurvLinear::BilinearForm, myNavierstokesTimeCurvLinear::LinearForm,
                      computeFreeSurfaceStress_onlyTP::LinearForm >
       timeSteppingProblem(V);

  // Create linear solver, set parameters and assign it to timeSteppingProblem
  dolfin::LinearSolver linear_solver;
  timeSteppingProblem.setLinearSolver (dolfin::reference_to_no_delete_pointer(linear_solver));


  Ivan::MovingTimeDependentProblem navierStokesProblem (dolfin::reference_to_no_delete_pointer (timeSteppingProblem),
                                                   t0,
                                                   dt, 
                                                   T, 
                                                   {"bilinear_form", "linear_form", "additional_form"},
                                                   {"bilinear_form", "linear_form"},
                                                   {"bilinear_form"}
                                                  );
         
std::cerr << "OK fino a " << __LINE__ << std::endl;

  // Boundary conditions
  dolfin::Constant noSlipDirichletBC (0.0, 0.0);
  dolfin::Constant freeSlipDirichletBC (0.0);
//inflow//  InflowDirichletData inflowDirichletBC;

  BottomBoundary inflowBoundary;
  LateralWall wallBoundary;
 // LateralInterior wallBoundaryInterior;
  TopBoundary freeSurface;
  TriplePointLeft triplePointLeft;
  TriplePointRight triplePointRight;
  TriplePointLeftVertex triplePointLeftVertex;
  TriplePointRightVertex triplePointRightVertex;

  //dolfin::FacetFunction<std::size_t> meshFacets (mesh);
  MyFacetFunction<std::size_t> meshFacets (mesh);
meshFacets.setPb (navierStokesProblem);
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
  navierStokesProblem.setIntegrationSubdomain ("bilinear_form",
        dolfin::reference_to_no_delete_pointer (meshFacets), dcp::SubdomainType::BOUNDARY_FACETS);
  navierStokesProblem.setIntegrationSubdomain ("linear_form",
        dolfin::reference_to_no_delete_pointer (meshFacets), dcp::SubdomainType::BOUNDARY_FACETS);
  navierStokesProblem.setIntegrationSubdomain ("additional_form",
        dolfin::reference_to_no_delete_pointer (additionalMeshFacets), dcp::SubdomainType::BOUNDARY_FACETS);
/*addVert//  navierStokesProblem.setIntegrationSubdomain ("additional_form",
        dolfin::reference_to_no_delete_pointer (additionalMeshVertices), dcp::SubdomainType::VERTICES);
*///addVert
std::cerr << "OK fino a " << __LINE__ << std::endl;
dolfin::plot(meshFacets, "meshFacets");dolfin::plot(additionalMeshFacets, "additionalMeshFacets");dolfin::interactive();

//inflow//  navierStokesProblem.addDirichletBC (
//inflow//        dolfin::DirichletBC (* (* navierStokesProblem.functionSpace())[0], inflowDirichletBC, inflowBoundary),
//inflow//        std::string("inflow"));
//  dcp::Subdomain inflowBoundary ((BottomBoundaryEvaluator ()));
//noInflow//  dcp::TimeDependentExpression inflowDirichletBC(2,(InflowDirichletBCEvaluator ()));
//noInflow//  navierStokesProblem.addTimeDependentDirichletBC (inflowDirichletBC, inflowBoundary, 0, "inflowBC");

std::unordered_map<std::size_t, double> bdval;
dolfin::DirichletBC dirBCwall (* (* (* navierStokesProblem.functionSpace())[0])[0], freeSlipDirichletBC, meshFacets,2, "geometric");
dirBCwall.get_boundary_values(bdval);
  navierStokesProblem.addDirichletBC (
          dirBCwall,
//        dolfin::DirichletBC (* (* (* navierStokesProblem.functionSpace())[0])[0], freeSlipDirichletBC, wallBoundary),
//        dolfin::DirichletBC (* (* navierStokesProblem.functionSpace())[0], noSlipDirichletBC, wallBoundary),
        std::string("wall"));
  #ifdef vanMourik
dolfin::DirichletBC dirBCin (* (* navierStokesProblem.functionSpace())[0], noSlipDirichletBC, meshFacets,3, "geometric");
  #else
dolfin::DirichletBC dirBCin (* (* (* navierStokesProblem.functionSpace())[0])[0], freeSlipDirichletBC, meshFacets,3, "geometric");
  #endif
dirBCin.get_boundary_values(bdval);
  navierStokesProblem.addDirichletBC (
          dirBCin,
//  #ifdef vanMourik
//            dolfin::DirichletBC (* (* navierStokesProblem.functionSpace())[0], noSlipDirichletBC, inflowBoundary),
//  #else
//        dolfin::DirichletBC (* (* (* navierStokesProblem.functionSpace())[0])[0], freeSlipDirichletBC, inflowBoundary),
//  #endif
        std::string("cut"));

  // Coefficients
  navierStokesProblem.setCoefficient ("bilinear_form", dolfin::reference_to_no_delete_pointer (nu),"nu");
  navierStokesProblem.setCoefficient ("bilinear_form", dolfin::reference_to_no_delete_pointer (beta),"beta");
  navierStokesProblem.setCoefficient ("linear_form", dolfin::reference_to_no_delete_pointer (gamma),"gamma");
  navierStokesProblem.setCoefficient ("linear_form", dolfin::reference_to_no_delete_pointer (stressAbove),"stressAbove");
  navierStokesProblem.setCoefficient ("linear_form", dolfin::reference_to_no_delete_pointer (stressBelow),"stressBelow");
  navierStokesProblem.setCoefficient ("linear_form", dolfin::reference_to_no_delete_pointer (*(new WallVelocity)),"wallVelocity");
  navierStokesProblem.setCoefficient ("linear_form", dolfin::reference_to_no_delete_pointer (gravityVersor),"gravityVersor");
  navierStokesProblem.setCoefficient ("linear_form", dolfin::reference_to_no_delete_pointer (beta),"beta");
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

