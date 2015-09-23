#include <dolfin.h>
#include <assert.h>
#include "MeshManager.h"
//#include "geometry.h"
#include "myNavierstokesTimeCurvLinear.h"
#include "UflToNewtonLinear.h"
#include "MovingTimeDependentProblem.h"
#include "utilities.h"
#include <dcp/subdomains/Subdomain.h>
#include "computeFreeSurfaceStress_onlyTP.h"

class TopBd : public dolfin::SubDomain
{
  public :
    bool inside (const dolfin::Array<double>& x, bool on_boundary) const
    {
        return on_boundary
               &&
               dolfin::near (x[1],ly);
    }
};
class BottomBd : public dolfin::SubDomain
{
  public :
    bool inside (const dolfin::Array<double>& x, bool on_boundary) const
    {
        return on_boundary
               &&
               dolfin::near (x[1],0);
    }
};
class LateralBd : public dolfin::SubDomain
{
  public :
    bool inside (const dolfin::Array<double>& x, bool on_boundary) const
    {
        return on_boundary
               &&
               ( dolfin::near (x[0],0) || dolfin::near (x[0],lx) );
    }
};

class InflowProfile : public dolfin::Expression
{
    void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x) const
    {
        values[0] = 0.0;
        values[1] = 0* 0.1 * 4/(1.1*lx*1.1*lx)*(x[0]+0.05*lx)*(1.05*lx-x[0]) ;
    }
    
    std::size_t value_rank() const
    {
      return 1;
    }

    std::size_t value_dimension(std::size_t i) const
    {
      return 2;
    }
};

int main (int argc, char * argv[])
{
 		dolfin::init(argc, argv);
dolfin::parameters["allow_extrapolation"] = true;

    // Mesh, labelling and mesh manager
 		dolfin::RectangleMesh mesh (0,0, lx,ly, nx,ny);
		dolfin::FacetFunction<std::size_t> meshFacets (dolfin::reference_to_no_delete_pointer (mesh));
	  meshFacets.set_all(0);
	  TopBd freeBoundary;
	  freeBoundary.mark (meshFacets, 1);
	  BottomBd inflowBoundary;
	  inflowBoundary.mark (meshFacets, 2);
	  LateralBd noSlipBoundary;
	  noSlipBoundary.mark (meshFacets, 3);
//	  LateralBd().mark (meshFacets, 3);
 		MeshManager<> meshManager (dolfin::reference_to_no_delete_pointer (mesh),
															dolfin::reference_to_no_delete_pointer (meshFacets), 0);
 		
 		// Function space, mesh manager setting and time stepping problem
 		myNavierstokesTimeCurvLinear::FunctionSpace V (* meshManager.mesh());
    meshManager.getDofOrder (V,0);
    TriplePointLeftVertex triplePointLeftVertex;
    TriplePointRightVertex triplePointRightVertex;
/*myNavierstokesTimeCurvLinear::BilinearForm form (V,V);
form.check();
    dolfin::VertexFunction<std::size_t> additionalMeshFunction (meshManager.mesh());
    additionalMeshFunction.set_all (0);
    triplePointLeftVertex.mark (additionalMeshFunction, 4);
    triplePointRightVertex.mark (additionalMeshFunction, 5);
for (dolfin::VertexIterator v (* meshManager.mesh()); !v.end(); ++v)
{
    if (triplePointLeftVertex.inside (dolfin::Array<double> (2,const_cast<double*>(v->x())),true))
      additionalMeshFunction[*v] = 4;
    if (triplePointRightVertex.inside (dolfin::Array<double> (2,const_cast<double*>(v->x())),true))
      additionalMeshFunction[*v] = 5;
}
//std::cerr << additionalMeshFunction.str(true) << std::endl;*/
const dolfin::GenericDofMap& dofmap (* V.dofmap());
// NB: qui sto usando la dofmap di un myNavierstokesTimeCurvLinear::FunctionSpace, mentre poi lo uso per processare un vettore che viene da computeFreeSurfaceStress_onlyTP::FunctionSpace
dolfin::la_index numDofs (dofmap.dofs().size());
const dolfin::GenericDofMap& subdofmap (* V[0]->dofmap());
std::vector<dolfin::la_index> velocityDofs (subdofmap.dofs());
std::vector<double> dofsCoords (dofmap.tabulate_all_coordinates (* meshManager.mesh()));
std::vector<dolfin::la_index> triplePointDofs;
assert (dofsCoords.size()==2*numDofs);
for (dolfin::la_index i=0; i!=numDofs; ++i)
{
  if ( triplePointLeftVertex.inside (dolfin::Array<double> (2, & dofsCoords[2*velocityDofs[i]]), true)
       ||
       triplePointRightVertex.inside (dolfin::Array<double> (2, & dofsCoords[2*velocityDofs[i]]), true) )
    triplePointDofs.push_back (velocityDofs[i]);
}
for (dolfin::la_index i=0; i!=triplePointDofs.size(); ++i) std::cerr << triplePointDofs[i] << ' '; std::cerr << std::endl;
for (dolfin::la_index i=0; i!=triplePointDofs.size(); ++i) std::cerr << dofsCoords[2*triplePointDofs[i]] << ' ' << dofsCoords[2*triplePointDofs[i]+1] << '\t'; std::cerr << std::endl;

 		Ivan::UflToNewtonLinear < myNavierstokesTimeCurvLinear::FunctionSpace, computeFreeSurfaceStress_onlyTP::FunctionSpace,
                      myNavierstokesTimeCurvLinear::BilinearForm, myNavierstokesTimeCurvLinear::LinearForm,
                      computeFreeSurfaceStress_onlyTP::LinearForm >
			  //timeSteppingProblem (dolfin::reference_to_no_delete_pointer (V));
/*			  timeSteppingProblem (dolfin::reference_to_no_delete_pointer(additionalMeshFunction),
                             {4,5},
                             dolfin::reference_to_no_delete_pointer (V));*/
        timeSteppingProblem (dolfin::reference_to_no_delete_pointer (V), triplePointDofs);
			  
    // Create linear solver, set parameters and assign it to timeSteppingProblem
		dolfin::LinearSolver linear_solver;
		timeSteppingProblem.setLinearSolver (dolfin::reference_to_no_delete_pointer(linear_solver));
		
		dolfin::Constant dt (0.00002);
	 	double t0 (0), T (t0+1000*dt);
//    meshManager.setImposedDisplacement (dolfin::reference_to_no_delete_pointer (V));

 		Ivan::MovingTimeDependentProblem navierStokesProblem (dolfin::reference_to_no_delete_pointer (meshManager),
		 																 										 dolfin::reference_to_no_delete_pointer (timeSteppingProblem),
                                                   			 t0,
                                                   			 dt, 
                                                   			 T, 
                                                   			 {"bilinear_form", "linear_form", "additional_form"},
                                                   			 {"bilinear_form", "linear_form"},
                                                   			 {"bilinear_form"}
                                                  			 );
    
    dolfin::plot (mesh,"out"); dolfin::interactive ();
    dolfin::plot (* meshManager.mesh(),"meshManager"); dolfin::interactive ();
/*    dolfin::plot (* V.mesh(),"V"); dolfin::interactive ();
    dolfin::plot (* timeSteppingProblem.mesh(),"timeSteppingProblem"); dolfin::interactive ();
    dolfin::plot (* timeSteppingProblem.functionSpace()->mesh(),"timeSteppingProblem.functionSpace"); dolfin::interactive ();
    dolfin::plot (* navierStokesProblem.mesh(),"NS pb"); dolfin::interactive ();
    dolfin::plot (* navierStokesProblem.functionSpace()->mesh(),"NS pb . functionSpace"); dolfin::interactive ();
    dolfin::plot (* navierStokesProblem.timeSteppingProblem().mesh(),"NS pb . timeSteppingProblem"); dolfin::interactive ();*/
    dolfin::plot (* navierStokesProblem.timeSteppingProblem().functionSpace()->mesh(),"NS pb . timeSteppingProblem.functionSpace"); dolfin::interactive ();

		// Setting the time-dependent problem
	  dolfin::Constant beta(0);
  	dolfin::Constant cosThetaS (cosThetaSVal);
    dolfin::Constant nu (2.081e-2 / rho);
		std::cerr << "Reynolds number = " << ustar*lx*rho/2.081e-2 << std::endl;
		dolfin::Constant gamma(1e-5);
  	dolfin::Constant stressBelow(0,0);
    InflowProfile inflowProfile;
 	  dcp::TimeDependentExpression stressAbove(2,(BoundaryStressEvaluator (0,0)));
    dolfin::Constant gravityVersor(0,-1);
    dolfin::Constant zeroVec (0,0);
    dolfin::Constant zero (0);
    dolfin::Constant t_partialOmega(0,1);
    navierStokesProblem.setCoefficient ("bilinear_form", dolfin::reference_to_no_delete_pointer (nu),"nu");
  	navierStokesProblem.setCoefficient ("bilinear_form", dolfin::reference_to_no_delete_pointer (beta),"beta");
  	navierStokesProblem.setCoefficient ("linear_form", dolfin::reference_to_no_delete_pointer (gamma),"gamma");
  	navierStokesProblem.setCoefficient ("linear_form", dolfin::reference_to_no_delete_pointer (stressAbove),"stressAbove");
  	navierStokesProblem.setCoefficient ("linear_form", dolfin::reference_to_no_delete_pointer (stressBelow),"stressBelow");
  	navierStokesProblem.setCoefficient ("linear_form", dolfin::reference_to_no_delete_pointer (zeroVec),"wallVelocity");
  	navierStokesProblem.setCoefficient ("linear_form", dolfin::reference_to_no_delete_pointer (gravityVersor),"gravityVersor");
  	navierStokesProblem.setCoefficient ("linear_form", dolfin::reference_to_no_delete_pointer (beta),"beta");
  	navierStokesProblem.setCoefficient ("additional_form", dolfin::reference_to_no_delete_pointer (gamma),"gamma");
  	navierStokesProblem.setCoefficient ("additional_form", dolfin::reference_to_no_delete_pointer (cosThetaS),"cosThetaS");
  	navierStokesProblem.setCoefficient ("additional_form", dolfin::reference_to_no_delete_pointer (t_partialOmega),"t_partialOmega");

		// Setting mesh labelling inside problem
		navierStokesProblem.setIntegrationSubdomain ("bilinear_form",
      dolfin::reference_to_no_delete_pointer (meshFacets), dcp::SubdomainType::BOUNDARY_FACETS);
		navierStokesProblem.setIntegrationSubdomain ("linear_form",
      dolfin::reference_to_no_delete_pointer (meshFacets), dcp::SubdomainType::BOUNDARY_FACETS);
		dolfin::FacetFunction<std::size_t> additionalMeshFacets (mesh);
		additionalMeshFacets.set_all (0);
		TriplePointLeft().mark(additionalMeshFacets,4);
		TriplePointRight().mark(additionalMeshFacets,5);
		navierStokesProblem.setIntegrationSubdomain ("additional_form",
      dolfin::reference_to_no_delete_pointer (additionalMeshFacets), dcp::SubdomainType::BOUNDARY_FACETS);

    // Setting Dirichlet BC
std::unordered_map<std::size_t, double> bdval;
//    dolfin::DirichletBC dirBCwall (* (* (* navierStokesProblem.functionSpace())[0])[0], zero, meshFacets, 3, "geometric");
   // dolfin::DirichletBC dirBCwall (* (* navierStokesProblem.functionSpace())[0], zeroVec, noSlipBoundary, "geometric");
    dolfin::DirichletBC dirBCwall (* (* (* navierStokesProblem.functionSpace())[0])[0], zero, noSlipBoundary, "geometric");
dirBCwall.get_boundary_values(bdval);
    navierStokesProblem.addDirichletBC (dirBCwall, std::string("wall"));
    dolfin::DirichletBC dirBCinflow (* (* navierStokesProblem.functionSpace())[0], inflowProfile, inflowBoundary, "geometric");
dirBCinflow.get_boundary_values(bdval);
    navierStokesProblem.addDirichletBC (dirBCinflow, std::string("inflow"));
    dolfin::DirichletBC dirBCoutflow (* (* (* navierStokesProblem.functionSpace())[0])[0], zero, freeBoundary, "geometric");
dirBCoutflow.get_boundary_values(bdval);
    navierStokesProblem.addDirichletBC (dirBCoutflow, std::string("outflow"));

/*XY xy;
dolfin::Function w (* meshManager.displacement()->function_space());
w = xy;
meshManager.moveMesh (w,"normal",1);

    dolfin::plot (mesh,"out"); dolfin::interactive ();
    dolfin::plot (* meshManager.mesh(),"meshManager"); dolfin::interactive ();
    dolfin::plot (* V.mesh(),"V"); dolfin::interactive ();
    dolfin::plot (* timeSteppingProblem.mesh(),"timeSteppingProblem"); dolfin::interactive ();
    dolfin::plot (* timeSteppingProblem.functionSpace()->mesh(),"timeSteppingProblem.functionSpace"); dolfin::interactive ();
    dolfin::plot (* navierStokesProblem.mesh(),"NS pb"); dolfin::interactive ();
    dolfin::plot (* navierStokesProblem.functionSpace()->mesh(),"NS pb . functionSpace"); dolfin::interactive ();
    dolfin::plot (* navierStokesProblem.timeSteppingProblem().mesh(),"NS pb . timeSteppingProblem"); dolfin::interactive ();
    dolfin::plot (* navierStokesProblem.timeSteppingProblem().functionSpace()->mesh(),"NS pb . timeSteppingProblem.functionSpace"); dolfin::interactive ();
    dolfin::plot (* (* navierStokesProblem.timeSteppingProblem().functionSpace())[0]->mesh(),"NS pb . timeSteppingProblem.functionSpace[0]"); dolfin::interactive ();
    dolfin::plot (* (* (* navierStokesProblem.timeSteppingProblem().functionSpace())[0])[1]->mesh(),"NS pb . timeSteppingProblem.functionSpace[0][1]"); dolfin::interactive ();

XYZ xyz;
dolfin::Function sol (V);
sol = xyz;*/

		navierStokesProblem.solve ();
dolfin::Function sol (navierStokesProblem.solution());
dolfin::plot (sol[0], "sol velocity");
dolfin::plot (sol[1], "sol pressure");
dolfin::plot (* sol.function_space()->mesh(), "final mesh"); dolfin::interactive();
//std::cerr << sol.vector()->str(true) << std::endl;

 		return 0;
}
