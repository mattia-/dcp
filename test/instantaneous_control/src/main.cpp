#include <dolfin.h>
#include <assert.h>
#include "navierstokes.h"
#include "navierstokesPreviousDomain.h"
#include "navierstokesAdjoint.h"
#include "navierstokesAdjointAdditional.h"
#include <dcp/differential_problems/MeshManager.h> //"MeshManager.h"
//#include "geometry.h"
#include <dcp/differential_problems/MovingLinearProblem.h>
#include <dcp/differential_problems/MovingTimeDependentProblem.h>
#include <dcp/differential_problems/utilities.h> //#include "utilities.h"
#include "InstantaneousControl.h" //<dcp/differential_problems/TimeDependentEquationSystem.h>
#include <dcp/subdomains/Subdomain.h>

#include "PostProcessor.h"
#include "computeFreeSurfaceStress_onlyTP.h"
#include "GetPot.h"
#include "laplaceALE.h"

#define EVITAPLOT
//#define SLOSHING
#define PARAB
//#define IMPOSEDISPL
  // TODO WARNING! Al momento funziona facendo setCoefficient di una Expression inflowProfile perche' inflowProfile non dipende da x[1]

//#define INITDISPL
//#define CLOSEDINFLOW

GetPot inputData;

ProblemData problemData;

class TPadditionalProcessor : public dcp::AdditionalProcessor
{
  public:

    void operator() (dolfin::Vector & processedVec, const dolfin::Vector & vec, const std::vector<dolfin::la_index> & additionalDofs) override
    {
      processedVec.setitem (additionalDofs.back(), (90==problemData.thetaS ? 0 : problemData.dt*problemData.gamma*cos(problemData.thetaS*3.14159265/180.0)*problemData.lx));
      
std::cerr << "processed values " << processedVec[additionalDofs.back()] << ", ";
    }
};

int main (int argc, char * argv[])
{
 		dolfin::init(argc, argv);
//dolfin::set_log_level (dolfin::DBG);
dolfin::parameters["allow_extrapolation"] = true;
    
    GetPot cmmLine (argc, argv);
    const std::string inputFileName = cmmLine.follow ("input.txt", 2, "-f", "--file");
    dolfin::info ("Reading input from file " + inputFileName);
    GetPot tmpGP (inputFileName.c_str());
    inputData = tmpGP;
    
    problemData.fileFrequency = inputData("fileFrequency", 10);
    problemData.saveDataInFile = inputData("saveDataInFile", false);
    problemData.savePvd = inputData("savePvd", false);
    problemData.startFromFlat = inputData("startFromFlat",true);
    problemData.inflowDirichlet = inputData("inflowDirichlet",false);

    problemData.lx = inputData("lx", 0.92e-3);
    problemData.ly = inputData("ly", problemData.lx);
    double yxratio (problemData.ly/problemData.lx);
    problemData.nx = inputData("nx", 16);
    problemData.ny = inputData("ny", 2*problemData.nx*yxratio);
		problemData.dt = inputData("dt", 0.0);
		dolfin::Constant dt (problemData.dt);
	 	double t0 (inputData("t0", 0.0)+(inputData("restartTimeStep",0)/2.0*problemData.fileFrequency*dt)),
           T  (inputData("T", t0));
    if (T <= t0)
    {
      std::size_t nT (inputData("nT", 0));
      T = t0 + nT*problemData.dt;
    }
std::cerr << "T  " << T << "   t0 " << t0;
      
    problemData.TP_EPS = problemData.lx/problemData.nx+6e-16;
    double  rho (inputData("rho", 1115.0)),
            re  (inputData("re", 0.0)),
            st  (inputData("st", 0.0)),
            ca  (inputData("ca", 0.0));
	  problemData.beta = inputData("betaSuH", 0.0);
	  dolfin::Constant beta (problemData.beta);
	  problemData.lambdastar = inputData("lambdastar", 0.5);
	  dolfin::Constant lambdastar (problemData.lambdastar);
	  problemData.noStokes = inputData("noStokes", 1.0);
	  dolfin::Constant noStokes (problemData.noStokes);
	  problemData.stabBulk = inputData("stabBulk", 1.0);
	  problemData.stabSigma = inputData("stabSigma", 0.0);
	  problemData.stabSGCL = inputData("stabSGCL", 0.0);
	  dolfin::Constant stabBulk (problemData.stabBulk);
	  dolfin::Constant stabSigma (problemData.stabSigma);
	  dolfin::Constant stabSGCL (problemData.stabSGCL);
	  dolfin::Constant stabPress (inputData("stabPress", 0.0));
  	problemData.thetaS = inputData("thetaS", 90.0);
  	dolfin::Constant cosThetaS (cos ( problemData.thetaS * 3.14159265/180.0 ));
    double mu (inputData("mu", 0.0));
    problemData.nu = mu / rho;
std::cerr << "mu  = " << mu << "     nu = " << problemData.nu << std::endl;
    dolfin::Constant nu (problemData.nu);
    double ustar (inputData("ustar", 0.1));
		std::cerr << "Reynolds number = " << ustar*problemData.lx*rho/mu << std::endl;
		problemData.gamma = inputData("gamma", 0.0);
		dolfin::Constant gamma(problemData.gamma);

//    const std::string saveRelPath = cmmLine.follow ("tmp/", 2, "-d", "--directory");
//    problemData.savepath = "/u/archive/dott/ifumagalli/dcp_test_output_continue/"+saveRelPath;
    problemData.savepath = cmmLine.follow ("tmp/", 2, "-d", "--directory");
    std::cerr << "Output directory: " + problemData.savepath << std::endl;

 		dolfin::RectangleMesh mesh (0,0, problemData.lx,problemData.ly, problemData.nx,problemData.ny);
		dolfin::FacetFunction<std::size_t> meshFacets (mesh);
	  meshFacets.set_all(0);
	  TopBd freeBoundary;
	  freeBoundary.mark (meshFacets, 1);
	  BottomBd inflowBoundary;
	  inflowBoundary.mark (meshFacets, 2);
/*	  LateralBd noSlipBoundary;
*/
    //AXISIMMETRIC
    RightBd noSlipBoundary;
	  noSlipBoundary.mark (meshFacets, 3);
    LeftBd slipBoundary;
    slipBoundary.mark (meshFacets, 4);
/* SIMPLE
//	  TopBd().mark (meshFacets, 1);
	  BottomBd fixedBoundary;
	  fixedBoundary.mark (meshFacets, 2);
//	  BottomBd().mark (meshFacets, 2);
	  LateralBd slipBoundary;
	  slipBoundary.mark (meshFacets, 3);
 */
//	  LateralBd().mark (meshFacets, 3);

    laplaceALE::FunctionSpace W (mesh);
    dcp::LinearProblem <laplaceALE::BilinearForm, laplaceALE::LinearForm> aleProblem (dolfin::reference_to_no_delete_pointer (W));
    dolfin::Constant uno (1.0), zero (0.0), penalty (1e30);
    aleProblem.setCoefficient ("bilinear_form",dolfin::reference_to_no_delete_pointer(uno),"k");
    //aleProblem.setCoefficient ("bilinear_form",dolfin::reference_to_no_delete_pointer(zero),"k");
    aleProblem.setCoefficient ("bilinear_form",dolfin::reference_to_no_delete_pointer(penalty),"penalty");
    aleProblem.setCoefficient ("linear_form",dolfin::reference_to_no_delete_pointer(uno),"k");
    //aleProblem.setCoefficient ("linear_form",dolfin::reference_to_no_delete_pointer(zero),"k");
    aleProblem.setCoefficient ("linear_form",dolfin::reference_to_no_delete_pointer(zero),"f");
    aleProblem.setCoefficient ("linear_form",dolfin::reference_to_no_delete_pointer(zero),"g");
    aleProblem.setCoefficient ("linear_form",dolfin::reference_to_no_delete_pointer(penalty),"penalty");
    
		aleProblem.setIntegrationSubdomain ("bilinear_form", dolfin::reference_to_no_delete_pointer(meshFacets), dcp::SubdomainType::BOUNDARY_FACETS);
		aleProblem.setIntegrationSubdomain ("linear_form", dolfin::reference_to_no_delete_pointer(meshFacets), dcp::SubdomainType::BOUNDARY_FACETS);
	  BottomBd bottomBoundary;
    dolfin::DirichletBC dirBCbottom (W, zero, bottomBoundary, "geometric");
std::unordered_map<std::size_t, double> bdval;
dirBCbottom.get_boundary_values(bdval);
    aleProblem.addDirichletBC (dirBCbottom, std::string("bottom"));
    laplaceALE::LinearForm::CoefficientSpace_u displCoeffSpace (mesh);
 		dcp::MeshManager<> meshManager (dolfin::reference_to_no_delete_pointer (mesh),
															dolfin::reference_to_no_delete_pointer (meshFacets),
                              dolfin::reference_to_no_delete_pointer (aleProblem),
                              displCoeffSpace,
                              0);

 		// Function space, mesh manager setting and time stepping problem
 		navierstokes::FunctionSpace V (* meshManager.mesh());
    meshManager.storeOrderedDofIdxs (V, {"ux","uy","p"});

//left    TriplePointLeftVertex triplePointLeftVertex;
    TriplePointRightVertex triplePointRightVertex;
/*navierstokes::BilinearForm form (V,V);
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
// NB: qui sto usando la dofmap di un navierstokes::FunctionSpace, mentre poi lo uso per processare un vettore che viene da computeFreeSurfaceStress_onlyTP::FunctionSpace
dolfin::la_index numDofs (dofmap.dofs().size());
const dolfin::GenericDofMap& subdofmap (* V[0]->dofmap());
std::vector<dolfin::la_index> velocityDofs (subdofmap.dofs());
std::vector<double> dofsCoords (dofmap.tabulate_all_coordinates (* meshManager.mesh()));
std::vector<dolfin::la_index> triplePointDofs;
assert (dofsCoords.size()==2*numDofs);
for (dolfin::la_index i=0; i!=numDofs; ++i)
{
  if ( //left  triplePointLeftVertex.inside (dolfin::Array<double> (2, & dofsCoords[2*velocityDofs[i]]), true)
       //left  ||
       triplePointRightVertex.inside (dolfin::Array<double> (2, & dofsCoords[2*velocityDofs[i]]), true) )
    triplePointDofs.push_back (velocityDofs[i]);
}
std::cerr << "triple points : "; for (dolfin::la_index i=0; i!=triplePointDofs.size(); ++i) std::cerr << triplePointDofs[i] << ' '; std::cerr << std::endl;
std::cerr << "                "; for (dolfin::la_index i=0; i!=triplePointDofs.size(); ++i) std::cerr << dofsCoords[2*triplePointDofs[i]] << ' ' << dofsCoords[2*triplePointDofs[i]+1] << '\t'; std::cerr << std::endl;
const dolfin::GenericDofMap& dofmap_w (* meshManager.displacement()->function_space()->dofmap());
dolfin::la_index numDofs_w (dofmap_w.dofs().size());
std::vector<dolfin::la_index> displacementDofs (dofmap_w.dofs());
std::vector<double> dofsCoords_w (dofmap_w.tabulate_all_coordinates (* meshManager.mesh()));
std::vector<dolfin::la_index> triplePointDofs_w;
assert (dofsCoords_w.size()==2*numDofs_w);
for (dolfin::la_index i=0; i!=numDofs_w; ++i)
{
  if ( //left  triplePointLeftVertex.inside (dolfin::Array<double> (2, & dofsCoords[2*velocityDofs[i]]), true)
       //left  ||
       triplePointRightVertex.inside (dolfin::Array<double> (2, & dofsCoords_w[2*displacementDofs[i]]), true) )
    triplePointDofs_w.push_back (displacementDofs[i]);
}
std::cerr << "w triple points : "; for (dolfin::la_index i=0; i!=triplePointDofs_w.size(); ++i) std::cerr << triplePointDofs_w[i] << ' '; std::cerr << std::endl;
std::cerr << "                "; for (dolfin::la_index i=0; i!=triplePointDofs_w.size(); ++i) std::cerr << dofsCoords_w[2*triplePointDofs_w[i]] << ' ' << dofsCoords_w[2*triplePointDofs_w[i]+1] << '\t'; std::cerr << std::endl;

    TPadditionalProcessor tpAdditionalProcessor;

 		dcp::MovingLinearProblem < navierstokes::FunctionSpace, computeFreeSurfaceStress_onlyTP::FunctionSpace,
                      navierstokes::BilinearForm, navierstokes::LinearForm,
                      navierstokesPreviousDomain::LinearForm, computeFreeSurfaceStress_onlyTP::LinearForm >
/*			  timeSteppingProblem (dolfin::reference_to_no_delete_pointer(additionalMeshFunction),
                             {4,5},
                             dolfin::reference_to_no_delete_pointer (V));*/
        timeSteppingProblem (dolfin::reference_to_no_delete_pointer (V), triplePointDofs, std::make_shared<TPadditionalProcessor> (tpAdditionalProcessor));
//REM//    timeSteppingProblem.setMeshManager (meshManager);
			  
    // Create linear solver, set parameters and assign it to timeSteppingProblem
		dolfin::LinearSolver linear_solver;
		timeSteppingProblem.setLinearSolver (dolfin::reference_to_no_delete_pointer(linear_solver));

#ifdef IMPOSEDISPL
    meshManager.setImposedDisplacement (dolfin::reference_to_no_delete_pointer (V));
#endif

#ifdef PARAB
 		dcp::MovingTimeDependentProblem navierStokesProblem (dolfin::reference_to_no_delete_pointer (meshManager),
		 																 										 dolfin::reference_to_no_delete_pointer (timeSteppingProblem),
                                                   			 t0,
                                                   			 problemData.dt, 
                                                   			 T, 
                                                   			 {"bilinear_form", "linear_form", "additional_form", "previous_form"},
                                                   			 {"bilinear_form", "linear_form", "previous_form"},
                                                   			 {"bilinear_form", "linear_form"}
                                                  			 );
#else
 		dcp::MovingTimeDependentProblem navierStokesProblem (dolfin::reference_to_no_delete_pointer (meshManager),
		 																 										 dolfin::reference_to_no_delete_pointer (timeSteppingProblem),
                                                   			 t0,
                                                   			 problemData.dt, 
                                                   			 T, 
                                                   			 {"bilinear_form", "linear_form", "additional_form", "previous_form"},
                                                   			 {"bilinear_form", "previous_form"},
                                                   			 {"bilinear_form"}
                                                  			 );
#endif

    navierStokesProblem.parameters["store_interval"] = 2;
    navierStokesProblem.parameters["time_stepping_solution_component"] = 0;
    
/*    dolfin::Constant initialSolution (0,0,0);
std::cerr << __LINE__ << std::endl;
    navierStokesProblem.setInitialSolution (initialSolution);
std::cerr << __LINE__ << std::endl;
    navierStokesProblem.setInitialSolution (initialSolution, 2);
*/

    //navierStokesProblem.setPostProcessor (new Ivan::PostProcessor(navierStokesProblem));

#ifndef EVITAPLOT
    dolfin::plot (mesh,"out"); dolfin::interactive ();
    dolfin::plot (* meshManager.mesh(),"meshManager"); dolfin::interactive ();
/*    dolfin::plot (* V.mesh(),"V"); dolfin::interactive ();
    dolfin::plot (* timeSteppingProblem.mesh(),"timeSteppingProblem"); dolfin::interactive ();
    dolfin::plot (* timeSteppingProblem.functionSpace()->mesh(),"timeSteppingProblem.functionSpace"); dolfin::interactive ();
    dolfin::plot (* navierStokesProblem.mesh(),"NS pb"); dolfin::interactive ();
    dolfin::plot (* navierStokesProblem.functionSpace()->mesh(),"NS pb . functionSpace"); dolfin::interactive ();
    dolfin::plot (* navierStokesProblem.timeSteppingProblem().mesh(),"NS pb . timeSteppingProblem"); dolfin::interactive ();*/
    dolfin::plot (* navierStokesProblem.timeSteppingProblem().functionSpace()->mesh(),"NS pb . timeSteppingProblem.functionSpace"); dolfin::interactive ();
#endif

		// Setting the time-dependent problem
  	dolfin::Constant stressBelow (inputData("zeta0_x" ,0.0), inputData("zeta0_y", 0.0));//9.81*ly);
 	  //DepthEvaluator<> stressBelowEvaluator (& meshManager);
 	  //dcp::TimeDependentExpression stressBelow(2, stressBelowEvaluator);
    InflowProfile inflowProfile;
    inflowProfile.setInflow (inputData("inflowVelocity",0.0));
 	  BoundaryStressEvaluator stressAboveEvaluator (inputData("stressAbove0_x", 0.0), - inputData("stressAbove0_y", 0.0));
 	  dcp::TimeDependentExpression stressAbove(2, stressAboveEvaluator);
    dolfin::Constant gravityVector(0,-inputData("gravityOn",1.0));
    dolfin::Constant zeroVec (0,0);
	  double wallVelY = inputData("wallVelY", 0.0);
    dolfin::Constant wallVel (0,wallVelY);
//ALREADY DECLARED//    dolfin::Constant zero (0);
    //dolfin::Constant t_partialOmega(0,1);
    dolfin::Constant t_partialOmega(0,problemData.lx); //cyl : integral_Gamma_w dGamma = 2*pi*R * integral_z dz
    //dolfin::Constant m_symmetryAxis(-1,0);
    navierStokesProblem.setCoefficient ("bilinear_form", dolfin::reference_to_no_delete_pointer (nu),"nu");
  	navierStokesProblem.setCoefficient ("bilinear_form", dolfin::reference_to_no_delete_pointer (beta),"beta");
  	navierStokesProblem.setCoefficient ("bilinear_form", dolfin::reference_to_no_delete_pointer (noStokes),"noStokes");
  	navierStokesProblem.setCoefficient ("bilinear_form", dolfin::reference_to_no_delete_pointer (stabBulk),"stabBulk");
  	navierStokesProblem.setCoefficient ("bilinear_form", dolfin::reference_to_no_delete_pointer (stabSigma),"stabSigma");
  	navierStokesProblem.setCoefficient ("bilinear_form", dolfin::reference_to_no_delete_pointer (stabSGCL),"stabSGCL");
  	navierStokesProblem.setCoefficient ("bilinear_form", dolfin::reference_to_no_delete_pointer (gamma),"gamma");
  	navierStokesProblem.setCoefficient ("bilinear_form", dolfin::reference_to_no_delete_pointer (stabPress),"stabPress");
#ifdef PARAB
    navierStokesProblem.setCoefficient ("linear_form", dolfin::reference_to_no_delete_pointer (nu),"nu");
#endif
  	navierStokesProblem.setCoefficient ("linear_form", dolfin::reference_to_no_delete_pointer (gamma),"gamma");
  	navierStokesProblem.setCoefficient ("linear_form", dolfin::reference_to_no_delete_pointer (stressAbove),"stressAbove");
  	navierStokesProblem.setCoefficient ("linear_form", dolfin::reference_to_no_delete_pointer (stressBelow),"stressBelow");
  	navierStokesProblem.setCoefficient ("linear_form", dolfin::reference_to_no_delete_pointer (wallVel),"wallVelocity");
  	navierStokesProblem.setCoefficient ("linear_form", dolfin::reference_to_no_delete_pointer (gravityVector),"gravityVector");
  	navierStokesProblem.setCoefficient ("linear_form", dolfin::reference_to_no_delete_pointer (beta),"beta");
  	navierStokesProblem.setCoefficient ("linear_form", dolfin::reference_to_no_delete_pointer (lambdastar),"lambdastar");
  	navierStokesProblem.setCoefficient ("linear_form", dolfin::reference_to_no_delete_pointer (stabBulk),"stabBulk");
  	//navierStokesProblem.setCoefficient ("linear_form", dolfin::reference_to_no_delete_pointer (stabPress),"stabPress");
  	navierStokesProblem.setCoefficient ("additional_form", dolfin::reference_to_no_delete_pointer (gamma),"gamma");
  	navierStokesProblem.setCoefficient ("additional_form", dolfin::reference_to_no_delete_pointer (cosThetaS),"cosThetaS");
  	navierStokesProblem.setCoefficient ("additional_form", dolfin::reference_to_no_delete_pointer (t_partialOmega),"t_partialOmega");
  	//navierStokesProblem.setCoefficient ("additional_form", dolfin::reference_to_no_delete_pointer (m_symmetryAxis),"m_symmetryAxis");
  	//navierStokesProblem.setCoefficient ("previous_form", dolfin::reference_to_no_delete_pointer (gamma),"gamma");
  	navierStokesProblem.setCoefficient ("previous_form", dolfin::reference_to_no_delete_pointer (stressAbove),"stressAbove");
  	navierStokesProblem.setCoefficient ("previous_form", dolfin::reference_to_no_delete_pointer (stressBelow),"stressBelow");
  	//navierStokesProblem.setCoefficient ("previous_form", dolfin::reference_to_no_delete_pointer (zeroVec),"wallVelocity");
  	navierStokesProblem.setCoefficient ("previous_form", dolfin::reference_to_no_delete_pointer (gravityVector),"gravityVector");
  	//navierStokesProblem.setCoefficient ("previous_form", dolfin::reference_to_no_delete_pointer (beta),"beta");
  	navierStokesProblem.setCoefficient ("previous_form", dolfin::reference_to_no_delete_pointer (lambdastar),"lambdastar");

/*    dolfin::Constant stabPress (inputData ("stabPress", 0.0) * problemData.lx*problemData.ly/(problemData.nx*problemData.ny));
          // it should be stabPress * h^2, but h=MaxCellEdgeLength is not available in the current FFC version
  	navierStokesProblem.setCoefficient ("bilinear_form", dolfin::reference_to_no_delete_pointer (stabPress),"stabPress");
  	navierStokesProblem.setCoefficient ("linear_form", dolfin::reference_to_no_delete_pointer (stabPress),"stabPress");*/
    dolfin::Constant rescale (inputData ("rescale", 1.0)); 
  	navierStokesProblem.setCoefficient ("linear_form", dolfin::reference_to_no_delete_pointer (rescale),"rescale");
#ifdef PARAB
    navierStokesProblem.setCoefficient ("bilinear_form", dolfin::reference_to_no_delete_pointer (inflowProfile), "uDir");
    navierStokesProblem.setCoefficient ("linear_form", dolfin::reference_to_no_delete_pointer (inflowProfile), "uDir");
    navierStokesProblem.setCoefficient ("previous_form", dolfin::reference_to_no_delete_pointer (inflowProfile), "uDir");
#endif

		// Setting mesh labelling inside problem
		navierStokesProblem.setIntegrationSubdomain ("bilinear_form",
      dolfin::reference_to_no_delete_pointer (meshFacets), dcp::SubdomainType::BOUNDARY_FACETS);
		navierStokesProblem.setIntegrationSubdomain ("linear_form",
      dolfin::reference_to_no_delete_pointer (meshFacets), dcp::SubdomainType::BOUNDARY_FACETS);
		dolfin::FacetFunction<std::size_t> additionalMeshFacets (mesh);
		additionalMeshFacets.set_all (0);
//left 	TriplePointLeft().mark(additionalMeshFacets,5);
		TriplePointRight().mark(additionalMeshFacets,6);
		navierStokesProblem.setIntegrationSubdomain ("additional_form",
      dolfin::reference_to_no_delete_pointer (additionalMeshFacets), dcp::SubdomainType::BOUNDARY_FACETS);

    // Setting Dirichlet BC
//ALREADY DECLARED// std::unordered_map<std::size_t, double> bdval;
//    dolfin::DirichletBC dirBCwall (* (* (* navierStokesProblem.functionSpace())[0])[0], zero, meshFacets, 3, "geometric");
   // dolfin::DirichletBC dirBCwall (* (* navierStokesProblem.functionSpace())[0], zeroVec, noSlipBoundary, "geometric");
    dolfin::DirichletBC dirBCwall (* (* (* navierStokesProblem.functionSpace())[0])[0], zero, noSlipBoundary, "geometric");
dirBCwall.get_boundary_values(bdval);
    navierStokesProblem.addDirichletBC (dirBCwall, std::string("wall"));
    dolfin::DirichletBC dirBCsymm (* (* (* navierStokesProblem.functionSpace())[0])[0], zero, slipBoundary, "geometric");
dirBCsymm.get_boundary_values(bdval);
    navierStokesProblem.addDirichletBC (dirBCsymm, std::string("symmetry"));
/*#ifdef SLOSHING
    dolfin::DirichletBC dirBCinflow (* (* navierStokesProblem.functionSpace())[0], zeroVec, inflowBoundary, "geometric");
dirBCinflow.get_boundary_values(bdval);
    navierStokesProblem.addDirichletBC (dirBCinflow, std::string("inflow"));
#else
  #ifdef PARAB
    dolfin::DirichletBC dirBCinflow (* (* navierStokesProblem.functionSpace())[0], zeroVec, inflowBoundary, "geometric");
  #else
    dolfin::DirichletBC dirBCinflow (* (* navierStokesProblem.functionSpace())[0], inflowProfile, inflowBoundary, "geometric");
  #endif
//    dolfin::DirichletBC dirBCinflow (* (* (* navierStokesProblem.functionSpace())[0])[0], zero, inflowBoundary, "geometric");
dirBCinflow.get_boundary_values(bdval);
    navierStokesProblem.addDirichletBC (dirBCinflow, std::string("inflow"));
#endif
*/
//#ifdef CLOSEDINFLOW
    if (problemData.inflowDirichlet)
    {
      dolfin::DirichletBC dirBCinflow (* (* navierStokesProblem.functionSpace())[0], zeroVec, inflowBoundary, "geometric");
dirBCinflow.get_boundary_values(bdval);
      navierStokesProblem.addDirichletBC (dirBCinflow, std::string("inflow"));
    }
//#endif

    dolfin::DirichletBC dirBCoutflow (* (* (* navierStokesProblem.functionSpace())[0])[0], zero, freeBoundary, "geometric");
dirBCoutflow.get_boundary_values(bdval);
    navierStokesProblem.addDirichletBC (dirBCoutflow, std::string("outflow"));

    // ------- adjoint problem ------ //

 		dcp::MovingLinearProblem < navierstokesAdjoint::FunctionSpace, navierstokesAdjointAdditional::FunctionSpace,
                      navierstokesAdjoint::BilinearForm, navierstokesAdjoint::LinearForm,
                      navierstokesAdjointAdditional::Form_previous, navierstokesAdjointAdditional::Form_additional >
        adjointProblem (dolfin::reference_to_no_delete_pointer (V));

    adjointProblem.setCoefficient ("bilinear_form", dolfin::reference_to_no_delete_pointer (dt),"dt");
    adjointProblem.setCoefficient ("bilinear_form", dolfin::reference_to_no_delete_pointer (nu),"nu");
  	adjointProblem.setCoefficient ("bilinear_form", dolfin::reference_to_no_delete_pointer (beta),"beta");
  	adjointProblem.setCoefficient ("bilinear_form", dolfin::reference_to_no_delete_pointer (noStokes),"noStokes");
  	adjointProblem.setCoefficient ("bilinear_form", dolfin::reference_to_no_delete_pointer (stabBulk),"stabBulk");
  	adjointProblem.setCoefficient ("bilinear_form", dolfin::reference_to_no_delete_pointer (stabSigma),"stabSigma");
  	adjointProblem.setCoefficient ("bilinear_form", dolfin::reference_to_no_delete_pointer (stabSGCL),"stabSGCL");
  	adjointProblem.setCoefficient ("bilinear_form", dolfin::reference_to_no_delete_pointer (gamma),"gamma");
  	adjointProblem.setCoefficient ("bilinear_form", dolfin::reference_to_no_delete_pointer (stabPress),"stabPress");
#ifdef PARAB
    adjointProblem.setCoefficient ("bilinear_form", dolfin::reference_to_no_delete_pointer (inflowProfile), "uDir");
#endif
    adjointProblem.setCoefficient ("linear_form", dolfin::reference_to_no_delete_pointer (dt),"dt");
    //?? MANCA u_old!!
//dolfin::Function u_old (* meshManager.displacement()); u_old = zeroVec; adjointProblem.setCoefficient ("linear_form",dolfin::reference_to_no_delete_pointer(u_old),"u_old");
//adjointProblem.setCoefficient ("bilinear_form",dolfin::reference_to_no_delete_pointer(zeroVec),"u_old");
std::cerr << __LINE__ << std::endl;
//    adjointProblem.setCoefficient ("linear_form",dolfin::reference_to_no_delete_pointer(zeroVec),"zeroVector");
    adjointProblem.setCoefficient ("additional_form",dolfin::reference_to_no_delete_pointer(zeroVec),"zeroVector");
    adjointProblem.setCoefficient ("previous_form",dolfin::reference_to_no_delete_pointer(zeroVec),"zeroVector");

		// Setting mesh labelling inside problem
		adjointProblem.setIntegrationSubdomain ("bilinear_form",
      dolfin::reference_to_no_delete_pointer (meshFacets), dcp::SubdomainType::BOUNDARY_FACETS);

    adjointProblem.addDirichletBC (dirBCwall, std::string("wall"));
    adjointProblem.addDirichletBC (dirBCsymm, std::string("symmetry"));
    adjointProblem.addDirichletBC (dirBCoutflow, std::string("outflow"));
    if (problemData.inflowDirichlet)
    {
      dolfin::DirichletBC dirBCinflow (* (* navierStokesProblem.functionSpace())[0], zeroVec, inflowBoundary, "geometric");
dirBCinflow.get_boundary_values(bdval);
      adjointProblem.addDirichletBC (dirBCinflow, std::string("inflow"));
    }

    // ------- system of equations ------- //

    navierstokesAdjointAdditional::Form_functional J (meshManager.mesh());
		J.set_exterior_facet_domains (dolfin::reference_to_no_delete_pointer (meshFacets));

    navierstokesAdjointAdditional::Form_integralNormal projector (meshManager.mesh());
    projector.set_exterior_facet_domains (dolfin::reference_to_no_delete_pointer (meshFacets));

  	dolfin::Constant control (0.0,0.0);
    dcp::InstantaneousControl tdSystem (dolfin::reference_to_no_delete_pointer (J), dolfin::reference_to_no_delete_pointer (projector), dolfin::reference_to_no_delete_pointer (control));
    tdSystem.parameters.add("control_name", "control");

    double penaltyFactor = inputData ("penaltyFactor", 1.0);
    tdSystem.parameters.add("penaltyFactor", penaltyFactor);

    double alpha = inputData ("alpha", 1.0);
    tdSystem.parameters.add("alpha", alpha);
    double alphaMin = inputData ("alphaMin", 1.0);
    tdSystem.parameters.add("alpha_min", alphaMin);

    tdSystem.setMeshManager (meshManager);

    tdSystem.addProblem ("primal", navierStokesProblem);
    tdSystem.addProblem ("adjoint", adjointProblem, false);
    tdSystem.addProblem ("ALE", aleProblem, false);

//     tdSystem.addLink ("adjoint", "primalSol", "linear_form", "primal");
//    tdSystem.addLinkToPreviousSolution ("adjoint", "u_old", "bilinear_form", "primal", 0, 1);
//    tdSystem.addLink ("adjoint", "w", "bilinear_form", "ALE");
//    tdSystem.addLink ("adjoint", "u", "linear_form", "primal", 0);
    tdSystem.addLink ("ALE", "u", "linear_form", "primal", 0);

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

// restart - mesh
    std::string restartMeshFilename (inputData ("restartMeshFile", ""));
    std::string restartTimeStep (inputData ("restartTimeStep", ""));

    if (! restartMeshFilename.empty())
    { 
      int rts (std::stoi(restartTimeStep));
      rts = (rts/2)*problemData.fileFrequency + (rts%2);
std::cerr << "restartTimeStep " << restartTimeStep;
      restartTimeStep = std::to_string (rts);
std::cerr << "       processed4hdf5 " << restartTimeStep << std::endl;

      restartMeshFilename += restartTimeStep; restartMeshFilename += ".hdf5";
      dolfin::Mesh restartMesh, restartMesh2;
      dolfin::HDF5File restartMeshFile (MPI_COMM_WORLD, restartMeshFilename, "r");
/*      restartMeshFile.read (restartMesh, std::to_string(5), true);
      restartMeshFile.read (restartMesh2, std::to_string(5), false);
//dolfin::plot (restartMesh, "restartMesh  5");
//dolfin::plot (restartMesh2, "restartMesh2  5");
      restartMeshFile.read (restartMesh, std::to_string(9), true);
      restartMeshFile.read (restartMesh2, std::to_string(9), false);
//dolfin::plot (restartMesh, "restartMesh  9");
//dolfin::plot (restartMesh2, "restartMesh2  9");
*/
      restartMeshFile.read (restartMesh, restartTimeStep, true);
//dolfin::plot (restartMesh, "mesh  restartTimeStep");
      restartMeshFile.read (restartMesh, restartTimeStep, false);
//dolfin::plot (restartMesh, "mesh2  restartTimeStep");
      restartMeshFile.close();
    
std::cerr << " OK fino a " << __LINE__ << std::endl;
      navierStokesProblem.meshManager().initializeMesh (restartMesh);
std::cerr << " OK fino a " << __LINE__ << std::endl;
      //dolfin::plot (* navierStokesProblem.mesh(), "initial problem mesh");
dolfin::interactive();
    }



    dolfin::Function tdsol (tdSystem.solution("primal"));
    tdSystem.solve();
    tdsol = tdSystem.solution("primal");
    std::cerr << std::endl << "Simulation completed!" << std::endl;
 		return 0;




dolfin::Function sol (navierStokesProblem.solution());
//dolfin::plot (sol[0], "sol velocity");
//dolfin::plot (sol[1], "sol pressure");
//dolfin::plot (* sol.function_space()->mesh(), "final mesh"); dolfin::interactive();
		navierStokesProblem.solve ();
#ifndef EVITAPLOT
//dolfin::Function sol (navierStokesProblem.solution());
sol = navierStokesProblem.solution();
dolfin::plot (sol[0], "sol velocity");
dolfin::plot (sol[1], "sol pressure");
dolfin::plot (* sol.function_space()->mesh(), "final mesh"); dolfin::interactive();
#endif
//std::cerr << sol.vector()->str(true) << std::endl;

    std::cerr << std::endl << "Simulation completed!" << std::endl;
 		return 0;
}
