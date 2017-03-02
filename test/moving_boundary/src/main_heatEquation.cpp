#include <dolfin.h>
#include <assert.h>
#include "myNavierstokesTimeCurvLinear.h"
#include "myNavierstokesTimeCurvLinearPreviousDomain.h"
#include <dcp/differential_problems/MeshManager.h>
#include <dcp/differential_problems/MovingLinearProblem.h>
#include "MovingTimeDependentProblem.h"
#include <dcp/differential_problems/utilities.h>
#include <dcp/subdomains/Subdomain.h>
#include "additionalForHeateq.h"
#include "GetPot.h"
#include "laplaceALEScal.h"
#include "heateq.h"
#include "heateqOld.h"

GetPot inputData;

ProblemData problemData;

class MyTopBd : public dolfin::SubDomain
{
  public :
    bool inside (const dolfin::Array<double>& x, bool on_boundary) const
    {
        return on_boundary
               &&
               dolfin::near (x[1],11);
    }
};
class MyBottomBd : public dolfin::SubDomain
{
  public :
    bool inside (const dolfin::Array<double>& x, bool on_boundary) const
    {
        return on_boundary
               &&
               dolfin::near (x[1],10);
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
    problemData.saveDataInFile = inputData("saveDataInFile", true);
    problemData.savePvd = inputData("savePvd", true);
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
      T = t0 + nT*dt;
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



  dolfin::RectangleMesh mesh (0,10, 1,11, problemData.nx,problemData.ny);
dolfin::plot(mesh, "initial mesh");
  laplaceALEScal::FunctionSpace W (mesh);
  dcp::LinearProblem <laplaceALEScal::BilinearForm, laplaceALEScal::LinearForm> aleProblem (dolfin::reference_to_no_delete_pointer (W));
  dolfin::Constant uno (1.0), zero (0.0), penalty (1e30);
    dolfin::Constant t_partialOmega(0,problemData.lx); //cyl : integral_Gamma_w dGamma = 2*pi*R * integral_z dz
  aleProblem.setCoefficient ("bilinear_form",dolfin::reference_to_no_delete_pointer(uno),"k");
  //aleProblem.setCoefficient ("bilinear_form",dolfin::reference_to_no_delete_pointer(zero),"k");
  aleProblem.setCoefficient ("bilinear_form",dolfin::reference_to_no_delete_pointer(penalty),"penalty");
  aleProblem.setCoefficient ("linear_form",dolfin::reference_to_no_delete_pointer(uno),"k");
  //aleProblem.setCoefficient ("linear_form",dolfin::reference_to_no_delete_pointer(zero),"k");
  aleProblem.setCoefficient ("linear_form",dolfin::reference_to_no_delete_pointer(zero),"f");
  aleProblem.setCoefficient ("linear_form",dolfin::reference_to_no_delete_pointer(zero),"g");
  aleProblem.setCoefficient ("linear_form",dolfin::reference_to_no_delete_pointer(penalty),"penalty");
    
  dolfin::FacetFunction<std::size_t> meshFacets (mesh);
	meshFacets.set_all(0);
	MyTopBd freeBoundary;
	freeBoundary.mark (meshFacets, 1);
	aleProblem.setIntegrationSubdomain ("bilinear_form", dolfin::reference_to_no_delete_pointer(meshFacets), dcp::SubdomainType::BOUNDARY_FACETS);
	aleProblem.setIntegrationSubdomain ("linear_form", dolfin::reference_to_no_delete_pointer(meshFacets), dcp::SubdomainType::BOUNDARY_FACETS);
	MyBottomBd bottomBoundary;
  dolfin::DirichletBC dirBCbottom (W, zero, bottomBoundary, "geometric");
std::unordered_map<std::size_t, double> bdval;
dirBCbottom.get_boundary_values(bdval);
  aleProblem.addDirichletBC (dirBCbottom, std::string("bottom"));

  laplaceALEScal::LinearForm::CoefficientSpace_u displCoeffSpace (mesh);

  dcp::MeshManager<> meshManager (dolfin::reference_to_no_delete_pointer (mesh),
 															    dolfin::reference_to_no_delete_pointer (meshFacets),
                                  dolfin::reference_to_no_delete_pointer (aleProblem),
                                  displCoeffSpace,
                                  0);
for (auto & e : meshManager.orderedDofs()) std::cerr << e.first << "   "; std::cerr << std::endl;
for (auto & e : meshManager.orderedDisplDofs()) std::cerr << e.first << "   "; std::cerr << std::endl;

  myNavierstokesTimeCurvLinear::FunctionSpace V (* meshManager.mesh());

/*  meshManager.getDofOrder (V, 0);
  std::vector<double> coords (V.dofmap()->tabulate_all_coordinates (* meshManager.mesh()));
{  const std::vector<std::size_t> & vec (meshManager.orderedDisplacementDofs ());
  std::cout << " ordDispl ";
  std::copy(vec.begin(), vec.end(), std::ostream_iterator<std::size_t>(std::cout, " "));
  std::cout << std::endl;
  std::cout << " coords ";
  std::copy(coords.begin(), coords.end(), std::ostream_iterator<double>(std::cout, " "));
  std::cout << std::endl;
}{  const std::vector<std::size_t> & vec (meshManager.orderedUXDofs ());
  std::cout << " ordUX ";
  std::copy(vec.begin(), vec.end(), std::ostream_iterator<std::size_t>(std::cout, " "));
  std::cout << std::endl;
}{  const std::vector<std::size_t> & vec (meshManager.orderedUYDofs ());
  std::cout << " ordUY ";
  std::copy(vec.begin(), vec.end(), std::ostream_iterator<std::size_t>(std::cout, " "));
  std::cout << std::endl;
}{  const std::vector<std::size_t> & vec (meshManager.orderedPressDofs ());
  std::cout << " ordP ";
  std::copy(vec.begin(), vec.end(), std::ostream_iterator<std::size_t>(std::cout, " "));
  std::cout << std::endl;
}{  const std::vector<std::size_t> & vec (meshManager.orderedWXDofs ());
  std::cout << " ordWX ";
  std::copy(vec.begin(), vec.end(), std::ostream_iterator<std::size_t>(std::cout, " "));
  std::cout << std::endl;
}{  const std::vector<std::size_t> & vec (meshManager.orderedWYDofs ());
  std::cout << " ordWY ";
  std::copy(vec.begin(), vec.end(), std::ostream_iterator<std::size_t>(std::cout, " "));
  std::cout << std::endl;
}
  const std::vector<std::pair<double, double> > & vecvec (meshManager.orderedDisplacementDofsCoords ());
  const std::vector<dolfin::la_index> & vec2 (meshManager.selectedOrderedUXDofs ());
  meshManager.selectedOrderedUYDofs ();
  meshManager.selectedOrderedPressDofs ();
  meshManager.selectedOrderedWXDofs ();
  meshManager.selectedOrderedWYDofs ();

  meshManager.storeOrderedDofIdxs (V, {"ux","uy","p"});

  dolfin::Function u (V);
  XYZ xyz;
  u = xyz;
  dcp::print2csv (u, problemData.savepath+"nome_",
           meshManager.orderedPressDofs().begin(),
           meshManager.orderedPressDofs().end(),
           u.function_space()->dofmap()->tabulate_all_coordinates (* meshManager.mesh())
           );
  meshManager.print2csv (u, problemData.savepath+"altro_", "A", true);

  dolfin::plot (u[0], "u");
  dolfin::plot (u[1], "p");
  dolfin::plot (* meshManager.displacement(), "w");
for (auto & e : meshManager.orderedDofs()) std::cerr << e.first << "   "; std::cerr << std::endl;
for (auto & e : meshManager.orderedDisplDofs()) std::cerr << e.first << "   "; std::cerr << std::endl;

dolfin::interactive();
*/

  heateq::FunctionSpace funSp (* meshManager.mesh());
  meshManager.storeOrderedDofIdxs (funSp, {""});
  dcp::MovingLinearProblem < heateq::FunctionSpace, additionalForHeateq::FunctionSpace,
                      heateq::BilinearForm, heateq::LinearForm,
                      heateqOld::LinearForm, additionalForHeateq::LinearForm >
        timeSteppingProblem (dolfin::reference_to_no_delete_pointer (funSp));
  
for (auto & e : meshManager.orderedDofs()) std::cerr << e.first << "   "; std::cerr << std::endl;
for (auto & e : meshManager.orderedDisplDofs()) std::cerr << e.first << "   "; std::cerr << std::endl;

	Ivan::MovingTimeDependentProblem heatEquation (dolfin::reference_to_no_delete_pointer (meshManager),
		 																 					  dolfin::reference_to_no_delete_pointer (timeSteppingProblem),
                                                0,
                                                problemData.dt,
                                                T, 
                                                {"bilinear_form", "linear_form", "additional_form"},
                                                {"previous_form"},
                                                {"bilinear_form"},
                                                -1
                                                );

  // UNNECESSARY : default value is -1
  //heatEquation.parameters["time_stepping_solution_component"] = -1; // "all" the components (i.e. the only one)

  heatEquation.setCoefficient ("bilinear_form", dolfin::reference_to_no_delete_pointer (uno),"nu");
  dolfin::Constant f (0.01);
  heatEquation.setCoefficient ("linear_form", dolfin::reference_to_no_delete_pointer (f),"f");
  heatEquation.setCoefficient ("linear_form", dolfin::reference_to_no_delete_pointer (zero),"g");
  heatEquation.setCoefficient ("additional_form", dolfin::reference_to_no_delete_pointer (zero),"gamma");

  dolfin::DirichletBC dirBC (* heatEquation.functionSpace(), zero, bottomBoundary, "geometric");
dirBC.get_boundary_values(bdval);
  heatEquation.addDirichletBC (dirBC, std::string("bottom"));

  heatEquation.solve ();

dolfin::interactive();

  return 0;
}
