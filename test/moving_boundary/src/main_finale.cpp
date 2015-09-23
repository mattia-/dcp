#include <dolfin.h>
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
               dolfin::near (x[1],2);
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
               ( dolfin::near (x[0],0) || dolfin::near (x[0],1) );
    }
};

int main (int argc, char * argv[])
{
 		dolfin::init(argc, argv);
dolfin::parameters["allow_extrapolation"] = true;

 		dolfin::RectangleMesh mesh (0,0, 1,2, 2,2);
		dolfin::FacetFunction<std::size_t> meshFacets (dolfin::reference_to_no_delete_pointer (mesh));
	  meshFacets.set_all(0);
/*	  TopBoundary freeBoundary;
	  freeBoundary.mark (meshFacets, 1);*/
	  TopBd().mark (meshFacets, 1);
/*	  BottomBoundary fixedBoundary;
	  fixedBoundary.mark (meshFacets, 2);*/
	  BottomBd().mark (meshFacets, 2);
/*	  LateralWall slipBoundary;
	  slipBoundary.mark (meshFacets, 3);*/
	  LateralBd().mark (meshFacets, 3);
 		MeshManager<> meshManager (dolfin::reference_to_no_delete_pointer (mesh),
															dolfin::reference_to_no_delete_pointer (meshFacets), 0);
 		
 		myNavierstokesTimeCurvLinear::FunctionSpace V (* meshManager.mesh());
 		Ivan::UflToNewtonLinear < myNavierstokesTimeCurvLinear::FunctionSpace, computeFreeSurfaceStress_onlyTP::FunctionSpace,
                      myNavierstokesTimeCurvLinear::BilinearForm, myNavierstokesTimeCurvLinear::LinearForm,
                      computeFreeSurfaceStress_onlyTP::LinearForm >
			  timeSteppingProblem (dolfin::reference_to_no_delete_pointer (V));
			  
    // Create linear solver, set parameters and assign it to timeSteppingProblem
		dolfin::LinearSolver linear_solver;
		timeSteppingProblem.setLinearSolver (dolfin::reference_to_no_delete_pointer(linear_solver));
		
		dolfin::Constant dt (0.00002);
	 	double t0 (0), T (t0+200*dt);
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

dolfin::Constant zero (1,2);
XY xy;
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
sol = xyz;
dolfin::plot (sol[0], "sol[0]");
std::cerr << sol.vector()->str(true) << std::endl;
std::vector<double> vals (sol.vector()->data(), sol.vector()->data()+sol.vector()->size());
for (auto it=vals.begin(); it!=vals.end(); it++) std::cerr << *it << ' '; std::cerr << std::endl;

 		return 0;
}
