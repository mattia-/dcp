//#include "MovingTimeDependentProblem.h"  // not required: it is already included by DefaultPostProcessor.h
#include "PostProcessor.h"
#include "discreteCurvature.h"
#include "discreteCurvature_onlyTP.h"
//addNotInMovingAbstract// #include "myNavierstokesTimeCurvLinear.h"
//addNotInMovingAbstract// #include "myNavierstokesTimeCurvLinearPreviousDomain.h"
#include "computeFreeSurfaceStress_onlyTP.h"
#include <dcp/differential_problems/utilities.h> //"utilities.h"

//#define GCL
  // if defined, assembles and file-writes the GCL dofs
  // NB: the initial files could be created anyway

//#define SINGLEBALANCEFILE

namespace Ivan
{

    PostProcessor::PostProcessor (dcp::MovingTimeDependentProblem & pb) :
      PostProcessor::PostProcessor (pb, problemData.savepath+"balanceTerms.csv", problemData.savepath+"intGCL.csv", problemData.savepath+"divGCL.csv", problemData.savepath+"selectedDofs.csv")
    {}

    PostProcessor::PostProcessor (dcp::MovingTimeDependentProblem & pb, std::string balanceFileName, std::string intGCLFileName, std::string divGCLFileName, std::string selectedDofsFileName) :
      DefaultPostProcessor (pb),
      formsDepOnSol_ (),
      formsDepOnOld_ (),
      formsDepOnDispl_ (),
      formsDepOnDir_ (),
      gclSpace1_ (* pb_.mesh()),
      gclSpace2_ (* pb_.mesh()),
      balanceFileName_ (balanceFileName),
      intGCLFileName_ (intGCLFileName),
      divGCLFileName_ (divGCLFileName),
      selectedDofsFileName_ (selectedDofsFileName),
/*      uxDofsNum_ ((2*problemData.nx+1)*(2*problemData.ny+1)), uyDofsNum_ (uxDofsNum_),
      pDofsNum_ ((problemData.nx+1)*(problemData.ny+1)),*/
      uxDofsNum_ (3*(3+5+3)), uyDofsNum_ (uxDofsNum_),
      pDofsNum_ (2*(2+3+2)),
      wxDofsNum_ (pDofsNum_), wyDofsNum_(pDofsNum_),
      uxDofsVals_ (new double [uxDofsNum_]), uyDofsVals_ (new double [uyDofsNum_]),
      pDofsVals_ (new double [pDofsNum_]),
      wxDofsVals_ (new double [wxDofsNum_]), wyDofsVals_ (new double [wyDofsNum_]),
      uxDofsIdxs_ (new dolfin::la_index [uxDofsNum_]), uyDofsIdxs_ (new dolfin::la_index [uyDofsNum_]),
      pDofsIdxs_ (new dolfin::la_index [pDofsNum_]),
      wxDofsIdxs_ (new dolfin::la_index [wxDofsNum_]), wyDofsIdxs_ (new dolfin::la_index [wyDofsNum_]),
      curvatureBilinearForm_ (new discreteCurvature::BilinearForm (* pb_.functionSpace(), * pb_.functionSpace())),
      curvatureLinearForm_ (new discreteCurvature::LinearForm (* pb_.functionSpace())),
      curvatureAdditionalForm_ (new discreteCurvature_onlyTP::LinearForm (* pb_.functionSpace())),
//addNotInMovingAbstract//      additionalFormDofs_ (std::static_pointer_cast<Ivan::MovingLinearProblem 
//addNotInMovingAbstract//                    < myNavierstokesTimeCurvLinear::FunctionSpace, computeFreeSurfaceStress_onlyTP::FunctionSpace,
//addNotInMovingAbstract//                      myNavierstokesTimeCurvLinear::BilinearForm, myNavierstokesTimeCurvLinear::LinearForm,
//addNotInMovingAbstract//                      myNavierstokesTimeCurvLinearPreviousDomain::LinearForm, computeFreeSurfaceStress_onlyTP::LinearForm > >
//addNotInMovingAbstract//                    (pb_.timeSteppingProblem_)->additionalFormDofs ()),
//addNotInMovingAbstract//      addIdxs_ (new dolfin::la_index [2]),
//addNotInMovingAbstract//      addIdxs_w_ (new dolfin::la_index [2]),
//addNotInMovingAbstract//      addVals_ (new double [2]),
//addNotInMovingAbstract//      addVals_w_ (new double [2]),
      auxiliary_ (pb_.functionSpace())
    {
      // coefficients storing
      std::shared_ptr<dolfin::GenericFunction> coeffPtr;
      coeffPtr.reset (new dolfin::Constant (1.0));
      coefficients_.emplace ("uno", coeffPtr);
      coeffPtr.reset (new dolfin::Constant (problemData.dt));
      coefficients_.emplace ("dt", coeffPtr);
      coeffPtr.reset (new dolfin::Constant (problemData.nu));
std::cerr << "checkNu " << problemData.nu << std::endl;
      coefficients_.emplace ("nu", coeffPtr);
      coeffPtr.reset (new dolfin::Constant (0,-inputData("gravityOn",1.0)));
      coefficients_.emplace ("gravityVector", coeffPtr);
  	  coeffPtr.reset (new dolfin::Constant (inputData("zeta0_x" ,0.0), inputData("zeta0_y", 0.0)));//9.81*ly);
      coefficients_.emplace ("stressBelow", coeffPtr);
      coeffPtr.reset (new dolfin::Constant (inputData("betaSuH",0.0)));
      coefficients_.emplace ("beta", coeffPtr);
      coeffPtr.reset (new dolfin::Constant (0,inputData("wallVelY",0.0)));
      coefficients_.emplace ("wallVelocity", coeffPtr);
      coeffPtr.reset (new dolfin::Constant (inputData("gamma",0.0)));
      coefficients_.emplace ("gamma", coeffPtr);
      coeffPtr.reset (new dolfin::Constant (cos(inputData("thetaS",90.0)*3.14159265/180.0)*problemData.lx));
      coefficients_.emplace ("cosThetaS", coeffPtr);
      coeffPtr.reset (new dolfin::Constant (inputData("stabPress",0.0)));
      coefficients_.emplace ("stabPress", coeffPtr);
      coeffPtr.reset (new dolfin::Constant (inputData("stabSGCL",0.0)));
      coefficients_.emplace ("stabSGCL", coeffPtr);

      // forms construction and storing
      std::shared_ptr<dolfin::Form> enPtr;

      enPtr.reset (new balanceTerms::Form_kinEn (* pb_.mesh()));
      formsNvalues_.emplace ("kinEn", std::make_pair (enPtr, double ()));
      formsNvalues_["kinEn"].first->set_coefficient ("dt", dolfin::reference_to_no_delete_pointer (* coefficients_["dt"]));
      formsDepOnSol_.push_back (formsNvalues_["kinEn"].first);

      enPtr.reset (new balanceTerms::Form_grav (* pb_.mesh()));
      formsNvalues_.emplace ("grav", std::make_pair (enPtr, double ()));
      formsNvalues_["grav"].first->set_coefficient ("dt", dolfin::reference_to_no_delete_pointer (* coefficients_["dt"]));
      formsNvalues_["grav"].first->set_coefficient ("f", dolfin::reference_to_no_delete_pointer (* coefficients_["gravityVector"])); 

      enPtr.reset (new balanceTerms::Form_gammaLength (* pb_.mesh()));
      formsNvalues_.emplace ("gammaLength", std::make_pair (enPtr, double ()));
      formsNvalues_["gammaLength"].first->set_coefficient ("dt", dolfin::reference_to_no_delete_pointer (* coefficients_["dt"])); 
      formsNvalues_["gammaLength"].first->set_coefficient ("gamma", dolfin::reference_to_no_delete_pointer (* coefficients_["gamma"])); 

      enPtr.reset (new balanceTerms::Form_sigmaLength (* pb_.mesh()));
      formsNvalues_.emplace ("sigmaLength", std::make_pair (enPtr, double ()));
      formsNvalues_["sigmaLength"].first->set_coefficient ("dt", dolfin::reference_to_no_delete_pointer (* coefficients_["dt"])); 
      formsNvalues_["sigmaLength"].first->set_coefficient ("gamma", dolfin::reference_to_no_delete_pointer (* coefficients_["gamma"])); 
      formsNvalues_["sigmaLength"].first->set_coefficient ("cosThetaS", dolfin::reference_to_no_delete_pointer (* coefficients_["cosThetaS"])); 

      enPtr.reset (new balanceTerms::Form_sigmaLength2 (* pb_.mesh()));
      formsNvalues_.emplace ("sigmaLength2", std::make_pair (enPtr, double ()));
      formsNvalues_["sigmaLength2"].first->set_coefficient ("uno", dolfin::reference_to_no_delete_pointer (* coefficients_["uno"])); 

      enPtr.reset (new balanceTerms::Form_kinEnFlux (* pb_.mesh()));
      formsNvalues_.emplace ("kinEnFlux", std::make_pair (enPtr, double ()));
      formsDepOnSol_.push_back (formsNvalues_["kinEnFlux"].first);
      formsDepOnOld_.push_back (formsNvalues_["kinEnFlux"].first);

      enPtr.reset (new balanceTerms::Form_gravFlux (* pb_.mesh()));
      formsNvalues_.emplace ("gravFlux", std::make_pair (enPtr, double ()));
      formsNvalues_["gravFlux"].first->set_coefficient ("f", dolfin::reference_to_no_delete_pointer (* coefficients_["gravityVector"])); 
      formsDepOnSol_.push_back (formsNvalues_["gravFlux"].first);

      enPtr.reset (new balanceTerms::Form_inflowStress (* pb_.mesh()));
      formsNvalues_.emplace ("inflowStress", std::make_pair (enPtr, double ()));
      formsNvalues_["inflowStress"].first->set_coefficient ("stressBelow", dolfin::reference_to_no_delete_pointer (* coefficients_["stressBelow"])); 
      formsDepOnSol_.push_back (formsNvalues_["inflowStress"].first);

      enPtr.reset (new balanceTerms::Form_viscPow (* pb_.mesh()));
      formsNvalues_.emplace ("viscPow", std::make_pair (enPtr, double ()));
      formsNvalues_["viscPow"].first->set_coefficient ("nu", dolfin::reference_to_no_delete_pointer (* coefficients_["nu"]));
      formsDepOnSol_.push_back (formsNvalues_["viscPow"].first);

      enPtr.reset (new balanceTerms::Form_nbc (* pb_.mesh()));
      formsNvalues_.emplace ("nbc", std::make_pair (enPtr, double ()));
      formsNvalues_["nbc"].first->set_coefficient ("beta", dolfin::reference_to_no_delete_pointer (* coefficients_["beta"])); 
      formsNvalues_["nbc"].first->set_coefficient ("wallVelocity", dolfin::reference_to_no_delete_pointer (* coefficients_["wallVelocity"])); 
      formsDepOnSol_.push_back (formsNvalues_["nbc"].first);

      enPtr.reset (new balanceTerms::Form_discr (* pb_.mesh()));
      formsNvalues_.emplace ("discr", std::make_pair (enPtr, double ()));
      formsNvalues_["discr"].first->set_coefficient ("dt", dolfin::reference_to_no_delete_pointer (* coefficients_["dt"]));
      formsDepOnSol_.push_back (formsNvalues_["discr"].first);
      formsDepOnOld_.push_back (formsNvalues_["discr"].first);

      enPtr.reset (new balanceTerms::Form_discrDiv (* pb_.mesh()));
      formsNvalues_.emplace ("discrDiv", std::make_pair (enPtr, double ()));
      formsNvalues_["discrDiv"].first->set_coefficient ("dt", dolfin::reference_to_no_delete_pointer (* coefficients_["dt"]));
      formsDepOnSol_.push_back (formsNvalues_["discrDiv"].first);
      formsDepOnOld_.push_back (formsNvalues_["discrDiv"].first);
      formsDepOnDispl_.push_back (formsNvalues_["discrDiv"].first);

      enPtr.reset (new balanceTerms::Form_epsg (* pb_.mesh()));
      formsNvaluesOnOld_.emplace ("epsg", std::make_pair (enPtr, double ()));
      formsNvaluesOnOld_["epsg"].first->set_coefficient ("dt", dolfin::reference_to_no_delete_pointer (* coefficients_["dt"]));
      formsNvaluesOnOld_["epsg"].first->set_coefficient ("f", dolfin::reference_to_no_delete_pointer (* coefficients_["gravityVector"])); 
      formsDepOnDispl_.push_back (formsNvaluesOnOld_["epsg"].first);

      enPtr.reset (new balanceTerms::Form_tgDivw (* pb_.mesh()));
      formsNvaluesOnOld_.emplace ("tgDivw", std::make_pair (enPtr, double ()));
      formsNvaluesOnOld_["tgDivw"].first->set_coefficient ("gamma", dolfin::reference_to_no_delete_pointer (* coefficients_["gamma"])); 
      formsNvaluesOnOld_["tgDivw"].first->set_coefficient ("dt", dolfin::reference_to_no_delete_pointer (* coefficients_["dt"]));
      formsDepOnDispl_.push_back (formsNvaluesOnOld_["tgDivw"].first);

      enPtr.reset (new balanceTerms::Form_tgDivSol (* pb_.mesh()));
      formsNvalues_.emplace ("tgDivSol", std::make_pair (enPtr, double ()));
      formsNvalues_["tgDivSol"].first->set_coefficient ("gamma", dolfin::reference_to_no_delete_pointer (* coefficients_["gamma"])); 
      formsDepOnSol_.push_back (formsNvalues_["tgDivSol"].first);

      enPtr.reset (new balanceTerms::Form_bdGravw (* pb_.mesh()));
      formsNvaluesOnOld_.emplace ("bdGravw", std::make_pair (enPtr, double ()));
      formsNvaluesOnOld_["bdGravw"].first->set_coefficient ("dt", dolfin::reference_to_no_delete_pointer (* coefficients_["dt"]));
      formsNvaluesOnOld_["bdGravw"].first->set_coefficient ("f", dolfin::reference_to_no_delete_pointer (* coefficients_["gravityVector"])); 
      formsDepOnDispl_.push_back (formsNvaluesOnOld_["bdGravw"].first);

      enPtr.reset (new balanceTerms::Form_bdGravSol (* pb_.mesh()));
      formsNvalues_.emplace ("bdGravSol", std::make_pair (enPtr, double ()));
      formsNvalues_["bdGravSol"].first->set_coefficient ("f", dolfin::reference_to_no_delete_pointer (* coefficients_["gravityVector"])); 
      formsDepOnSol_.push_back (formsNvalues_["bdGravSol"].first);

      enPtr.reset (new balanceTerms::Form_supg (* pb_.mesh()));
      formsNvalues_.emplace ("supg", std::make_pair (enPtr, double ()));
      formsNvalues_["supg"].first->set_coefficient ("stabPress", dolfin::reference_to_no_delete_pointer (* coefficients_["stabPress"])); 
      formsDepOnSol_.push_back (formsNvalues_["supg"].first);

      enPtr.reset (new balanceTerms::Form_gravPow (* pb_.mesh()));
      formsNvalues_.emplace ("gravPow", std::make_pair (enPtr, double ()));
      formsNvalues_["gravPow"].first->set_coefficient ("f", dolfin::reference_to_no_delete_pointer (* coefficients_["gravityVector"])); 
      formsDepOnSol_.push_back (formsNvalues_["gravPow"].first);

      enPtr.reset (new balanceTerms::Form_gDivSol (* pb_.mesh()));
      formsNvalues_.emplace ("gDivSol", std::make_pair (enPtr, double ()));
      formsNvalues_["gDivSol"].first->set_coefficient ("f", dolfin::reference_to_no_delete_pointer (* coefficients_["gravityVector"])); 
      formsDepOnSol_.push_back (formsNvalues_["gDivSol"].first);

      enPtr.reset (new balanceTerms::Form_tgDivSigmaw (* pb_.mesh()));
      formsNvaluesOnOld_.emplace ("tgDivSigmaw", std::make_pair (enPtr, double ()));
      formsNvaluesOnOld_["tgDivSigmaw"].first->set_coefficient ("dt", dolfin::reference_to_no_delete_pointer (* coefficients_["dt"])); 
      formsNvaluesOnOld_["tgDivSigmaw"].first->set_coefficient ("gamma", dolfin::reference_to_no_delete_pointer (* coefficients_["gamma"])); 
      formsNvaluesOnOld_["tgDivSigmaw"].first->set_coefficient ("cosThetaS", dolfin::reference_to_no_delete_pointer (* coefficients_["cosThetaS"])); 
      formsDepOnDispl_.push_back (formsNvaluesOnOld_["tgDivSigmaw"].first);

      enPtr.reset (new balanceTerms::Form_SGCLstab (* pb_.mesh()));
      formsNvaluesOnOld_.emplace ("SGCLstab", std::make_pair (enPtr, double ()));
      formsNvaluesOnOld_["SGCLstab"].first->set_coefficient ("stabSGCL", dolfin::reference_to_no_delete_pointer (* coefficients_["stabSGCL"])); 
      formsNvaluesOnOld_["SGCLstab"].first->set_coefficient ("gamma", dolfin::reference_to_no_delete_pointer (* coefficients_["gamma"])); 
      formsNvaluesOnOld_["SGCLstab"].first->set_coefficient ("dt", dolfin::reference_to_no_delete_pointer (* coefficients_["dt"])); 
      formsDepOnDispl_.push_back (formsNvaluesOnOld_["SGCLstab"].first);

      enPtr.reset (new balanceTerms::Form_SGCLstabOrig (* pb_.mesh()));
      formsNvalues_.emplace ("SGCLstabOrig", std::make_pair (enPtr, double ()));
std::cerr << " OK fino a " << __LINE__ << std::endl;
      formsNvalues_["SGCLstabOrig"].first->set_coefficient ("stabSGCL", dolfin::reference_to_no_delete_pointer (* coefficients_["stabSGCL"])); 
std::cerr << " OK fino a " << __LINE__ << std::endl;
      formsNvalues_["SGCLstabOrig"].first->set_coefficient ("gamma", dolfin::reference_to_no_delete_pointer (* coefficients_["gamma"])); 
std::cerr << " OK fino a " << __LINE__ << std::endl;
      //formsNvalues_["SGCLstabOrig"].first->set_coefficient ("dt", dolfin::reference_to_no_delete_pointer (* coefficients_["dt"])); 
std::cerr << " OK fino a " << __LINE__ << std::endl;
      formsDepOnSol_.push_back (formsNvalues_["SGCLstabOrig"].first);
      formsDepOnDispl_.push_back (formsNvalues_["SGCLstabOrig"].first);

#ifdef GCL
      enPtr.reset (new balanceTerms::Form_intGCL1 (gclSpace1_));
      gclForms_.emplace ("intGCL1", enPtr);

      enPtr.reset (new balanceTerms::Form_divGCL1 (gclSpace1_));
      gclForms_.emplace ("divGCL1", enPtr);
      formsDepOnDispl_.push_back (gclForms_["divGCL1"]);

      enPtr.reset (new balanceTerms::Form_intGCL2 (gclSpace2_));
      gclForms_.emplace ("intGCL2", enPtr);

      enPtr.reset (new balanceTerms::Form_divGCL2 (gclSpace2_));
      gclForms_.emplace ("divGCL2", enPtr);
      formsDepOnDispl_.push_back (gclForms_["divGCL2"]);
#endif

      // Dirichlet datum and meshFunction setting
      std::shared_ptr<const dolfin::MeshFunction<std::size_t> > meshFunction (pb_.meshManager().meshFacets());
      for (auto & formNvalue : formsNvalues_)
        formNvalue.second.first -> set_exterior_facet_domains (meshFunction);
      for (auto & formNvalue : formsNvaluesOnOld_)
        formNvalue.second.first -> set_exterior_facet_domains (meshFunction);

#ifdef SINGLEBALANCEFILE
      std::ofstream balanceFile (balanceFileName_, std::ofstream::out);
      balanceFile << "kinEn, grav, gammaLength, sigmaLength, kinEnFlux, gravFlux, inflowStress, viscPow, nbc, discr, discrDiv, epsg, tgDivw, tgDivSol, bdGravw, bdGravSol, supg, tpW, tpVelo, gravPow, gDivSol, tgDivSigmaw, SGCLstab, SGCLstabOrig";
      balanceFile << std::endl;
      balanceFile.close();
#endif

      uxDofsIdxs_ = pb_.meshManager().selectedOrderedUXDofs().data(); uyDofsIdxs_ = pb_.meshManager().selectedOrderedUYDofs().data();
      pDofsIdxs_ = pb_.meshManager().selectedOrderedPressDofs().data();
      wxDofsIdxs_ = pb_.meshManager().selectedOrderedWXDofs().data(); wyDofsIdxs_ = pb_.meshManager().selectedOrderedWYDofs().data();

const dolfin::GenericDofMap& dofmap_w (* pb_.meshManager().displacement()->function_space()->dofmap());
dolfin::la_index numDofs_w (dofmap_w.dofs().size());
std::vector<dolfin::la_index> displacementDofs (dofmap_w.dofs());
std::vector<double> dofsCoords_w (dofmap_w.tabulate_all_coordinates (* pb_.meshManager().mesh()));
std::vector<dolfin::la_index> triplePointDofs_w;
assert (dofsCoords_w.size()==2*numDofs_w);
TriplePointRightVertex triplePointRightVertex;
for (dolfin::la_index i=0; i!=numDofs_w; ++i)
{
  if ( //left  triplePointLeftVertex.inside (dolfin::Array<double> (2, & dofsCoords[2*velocityDofs[i]]), true)
       //left  ||
       triplePointRightVertex.inside (dolfin::Array<double> (2, & dofsCoords_w[2*displacementDofs[i]]), true) )
    triplePointDofs_w.push_back (displacementDofs[i]);
}
//addNotInMovingAbstract//      for (std::size_t i (0); i<2; ++i)
//addNotInMovingAbstract//      {
//addNotInMovingAbstract//        addIdxs_[i] = additionalFormDofs_[i];
//addNotInMovingAbstract//        addIdxs_w_[i] = triplePointDofs_w[i];
//addNotInMovingAbstract//      }
//addNotInMovingAbstract//std::cerr << "PostProcessor::addIdxs_:addIdxs_w_ "; for (std::size_t iii(0); iii<2; ++iii) std::cerr << addIdxs_[iii] << ':' << addIdxs_w_[iii] << ' '; std::cerr << std::endl;
    }


    void PostProcessor::onOldDomain (int timeStep, const dcp::MovingTimeDependentProblem * const pb)
    {
      std::vector<double> bCoords1 (gclSpace1_.dofmap()->tabulate_all_coordinates (* pb_.mesh()));
      std::vector<std::size_t> bIdxs1 (bCoords1.size()/2);
      std::iota (bIdxs1.begin(), bIdxs1.end(), 0);
      dolfin::Vector bdiv1;
      std::vector<double> bCoords2 (gclSpace2_.dofmap()->tabulate_all_coordinates (* pb_.mesh()));
      std::vector<std::size_t> bIdxs2 (bCoords2.size()/2);
      std::iota (bIdxs2.begin(), bIdxs2.end(), 0);
      dolfin::Vector bdiv2;

      const dolfin::Function & sol (pb_.solution()),
                          & oldSol (pb_.solutions().size()>1 ? pb_.solutions()[pb_.solutions().size()-2].second : sol),
                               & w (* pb_.meshManager().displacement());
      for (auto & form : formsDepOnSol_)
        //form -> set_coefficient ("u", dolfin::reference_to_no_delete_pointer (sol));
        form -> set_coefficient ("sol", dolfin::reference_to_no_delete_pointer (sol));
      for (auto & form : formsDepOnOld_)
        //form -> set_coefficient ("u_old", dolfin::reference_to_no_delete_pointer (oldSol));
        form -> set_coefficient ("sol_old", dolfin::reference_to_no_delete_pointer (oldSol));
      for (auto & form : formsDepOnDispl_)
        form -> set_coefficient ("w", dolfin::reference_to_no_delete_pointer (w));

      for (auto & formNvalue : formsNvaluesOnOld_)
      {
//        std::cerr << "processing " << formNvalue.first << std::endl;
        formNvalue.second.second = dolfin::assemble (* formNvalue.second.first);
      }
#ifdef GCL
      dolfin::assemble (bdiv1, * gclForms_["divGCL1"]);
      dolfin::Function bdivFun1 (gclSpace1_);
      * bdivFun1.vector() = bdiv1;
      dcp::print2csv (bdivFun1, divGCLFileName_+"old1."+std::to_string(timeStep), bIdxs1.begin(), bIdxs1.end(), bCoords1);

      dolfin::assemble (bdiv2, * gclForms_["divGCL2"]);
      dolfin::Function bdivFun2 (gclSpace2_);
      * bdivFun2.vector() = bdiv2;
      dcp::print2csv (bdivFun2, divGCLFileName_+"old2."+std::to_string(timeStep), bIdxs2.begin(), bIdxs2.end(), bCoords2);
#endif
    }

    void PostProcessor::operator() (int timeStep, const dcp::MovingTimeDependentProblem * const pb)
    {
      const dolfin::Function & sol (pb_.solution()),
                          & oldSol (pb_.solutions().size()>1 ? pb_.solutions()[pb_.solutions().size()-2].second : sol),
                               & w (* pb_.meshManager().displacement());
      for (auto & form : formsDepOnSol_)
        //form -> set_coefficient ("u", dolfin::reference_to_no_delete_pointer (sol));
        form -> set_coefficient ("sol", dolfin::reference_to_no_delete_pointer (sol));
      for (auto & form : formsDepOnOld_)
        //form -> set_coefficient ("u_old", dolfin::reference_to_no_delete_pointer (oldSol));
        form -> set_coefficient ("sol_old", dolfin::reference_to_no_delete_pointer (oldSol));
      for (auto & form : formsDepOnDispl_)
        form -> set_coefficient ("w", dolfin::reference_to_no_delete_pointer (w));

      for (auto & formNvalue : formsNvalues_)
      {
        std::cerr << "processing " << formNvalue.first << std::endl;
        formNvalue.second.second = dolfin::assemble (* formNvalue.second.first);
      }

#ifdef SINGLEBALANCEFILE
      std::ofstream balanceFile (balanceFileName_, std::ofstream::app);
#else
      std::ofstream balanceFile (balanceFileName_+"."+std::to_string(timeStep), std::ofstream::out);
#endif
      balanceFile.setf ( std::ios::scientific, std::ios::floatfield);
      balanceFile.precision (15);
//      balanceFile << "kinEn, grav, gammaLength, sigmaLength, kinEnFlux, gravFlux, inflowStress, viscPow, nbc, discr, discrDiv, epsg, tgDivw, tgDivSol, bdGravw, bdGravSol, supg, tpW, tpVelo, gravPow, gDivSol, tgDivSigmaw, SGCLstab, SGCLstabOrig";
      balanceFile << formsNvalues_["kinEn"].second << ',';
      balanceFile << formsNvalues_["grav"].second << ',';
      balanceFile << formsNvalues_["gammaLength"].second << ',';
      balanceFile << formsNvalues_["sigmaLength"].second << ',';
      balanceFile << formsNvalues_["kinEnFlux"].second << ',';
      balanceFile << formsNvalues_["gravFlux"].second << ',';
      balanceFile << formsNvalues_["inflowStress"].second << ',';
      balanceFile << formsNvalues_["viscPow"].second << ',';
      balanceFile << formsNvalues_["nbc"].second << ',';
      balanceFile << formsNvalues_["discr"].second << ',';
      balanceFile << formsNvalues_["discrDiv"].second << ',';
      balanceFile << formsNvaluesOnOld_["epsg"].second << ',';
      //outputFile << formsNvalues_["epsg"].second << ',';
      balanceFile << formsNvaluesOnOld_["tgDivw"].second << ',';
      balanceFile << formsNvalues_["tgDivSol"].second << ',';
      balanceFile << formsNvaluesOnOld_["bdGravw"].second << ',';
      balanceFile << formsNvalues_["bdGravSol"].second << ',';
      balanceFile << formsNvalues_["supg"].second << ',';
std::vector<double> blabladof (w.function_space()->dofmap()->tabulate_all_coordinates(* w.function_space()->mesh()));
//addNotInMovingAbstract//std::cerr << "DOF w " << blabladof[2*addIdxs_w_[0]] << ',' << blabladof[2*addIdxs_w_[0]+1] << ' ' << blabladof[2*addIdxs_w_[1]] << ',' << blabladof[2*addIdxs_w_[1]+1] << " -> " << (* w.vector())[addIdxs_w_[1]] << std::endl;
//addNotInMovingAbstract//      balanceFile << problemData.gamma*cos(problemData.thetaS*3.14159265/180.0)*problemData.lx*(* w.vector())[addIdxs_w_[1]]/problemData.dt << ',';
//addNotInMovingAbstract//      balanceFile << problemData.gamma*cos(problemData.thetaS*3.14159265/180.0)*problemData.lx*(* sol.vector())[addIdxs_[1]] << ',';
      balanceFile << formsNvalues_["gravPow"].second << ',';
      balanceFile << formsNvalues_["gDivSol"].second << ',';
      balanceFile << formsNvaluesOnOld_["tgDivSigmaw"].second << ',';
      balanceFile << formsNvaluesOnOld_["SGCLstab"].second << ',';
      balanceFile << formsNvalues_["SGCLstabOrig"].second << ',';
      balanceFile << formsNvalues_["sigmaLength2"].second << std::endl;
      balanceFile.close();

      std::vector<double> bCoords1 (gclSpace1_.dofmap()->tabulate_all_coordinates (* pb_.mesh()));
      std::vector<std::size_t> bIdxs1 (bCoords1.size()/2);
      std::iota (bIdxs1.begin(), bIdxs1.end(), 0);
      std::vector<double> bCoords2 (gclSpace2_.dofmap()->tabulate_all_coordinates (* pb_.mesh()));
      std::vector<std::size_t> bIdxs2 (bCoords2.size()/2);
      std::iota (bIdxs2.begin(), bIdxs2.end(), 0);
#ifdef GCL
      dolfin::Vector b1;
      dolfin::assemble (b1, * gclForms_["intGCL1"]);
      dolfin::Function bFun1 (gclSpace1_);
      * bFun1.vector() = b1;
      dcp::print2csv (bFun1, intGCLFileName_+"1."+std::to_string(timeStep), bIdxs1.begin(), bIdxs1.end(), bCoords1);
      dolfin::Vector bdiv1;
      dolfin::assemble (bdiv1, * gclForms_["divGCL1"]);
      * bFun1.vector() = bdiv1;
      dcp::print2csv (bFun1, divGCLFileName_+"1."+std::to_string(timeStep), bIdxs1.begin(), bIdxs1.end(), bCoords1);
      dolfin::Vector b2;
      dolfin::assemble (b2, * gclForms_["intGCL2"]);
      dolfin::Function bFun2 (gclSpace2_);
      * bFun2.vector() = b2;
      dcp::print2csv (bFun2, intGCLFileName_+"2."+std::to_string(timeStep), bIdxs2.begin(), bIdxs2.end(), bCoords2);
      dolfin::Vector bdiv2;
      dolfin::assemble (bdiv2, * gclForms_["divGCL2"]);
      * bFun2.vector() = bdiv2;
      dcp::print2csv (bFun2, divGCLFileName_+"2."+std::to_string(timeStep), bIdxs2.begin(), bIdxs2.end(), bCoords2);
#endif

//std::cerr<<"uxDofsIdxs ("<<uxDofsNum_<< ") "; for (dolfin::la_index i (0); i<uxDofsNum_; ++i) std::cerr << uxDofsIdxs_[i] << ", "; std::cerr<<std::endl;
//std::cerr<<"uyDofsIdxs ("<<uyDofsNum_<< ") "; for (dolfin::la_index i (0); i<uyDofsNum_; ++i) std::cerr << uyDofsIdxs_[i] << ", "; std::cerr<<std::endl;
//std::cerr<<"pDofsIdxs ("<<pDofsNum_<< ") "; for (dolfin::la_index i (0); i<pDofsNum_; ++i) std::cerr << pDofsIdxs_[i] << ", "; std::cerr<<std::endl;
//std::cerr<<"wxDofsIdxs ("<<wxDofsNum_<< ") "; for (dolfin::la_index i (0); i<wxDofsNum_; ++i) std::cerr << wxDofsIdxs_[i] << ", "; std::cerr<<std::endl;
//std::cerr<<"wyDofsIdxs ("<<wyDofsNum_<< ") "; for (dolfin::la_index i (0); i<wyDofsNum_; ++i) std::cerr << wyDofsIdxs_[i] << ", "; std::cerr<<std::endl;

//limitOutputFiles//      sol.vector()->get (uxDofsVals_, & uxDofsNum_, & uxDofsIdxs_);
//limitOutputFiles//      sol.vector()->get (uyDofsVals_, & uyDofsNum_, & uyDofsIdxs_);
//limitOutputFiles//      sol.vector()->get (pDofsVals_, & pDofsNum_, & pDofsIdxs_);
//limitOutputFiles//      w.vector()->get (wxDofsVals_, & wxDofsNum_, & wxDofsIdxs_);
//limitOutputFiles//      w.vector()->get (wyDofsVals_, & wyDofsNum_, & wyDofsIdxs_);
//limitOutputFiles//      print2csv (std::vector<const double *> ({uxDofsVals_, uyDofsVals_, pDofsVals_, wxDofsVals_, wyDofsVals_}),
//limitOutputFiles//                 std::vector<const dolfin::la_index *> ({& uxDofsNum_, & uyDofsNum_, & pDofsNum_, & wxDofsNum_, & wyDofsNum_}),
//limitOutputFiles//                 selectedDofsFileName_+"."+std::to_string(timeStep));
//limitOutputFiles//
//limitOutputFiles//      //discrete curvature
//limitOutputFiles//
//limitOutputFiles//      dolfin::Assembler assembler;
//limitOutputFiles//      assembler.assemble (this->curvatureMatrix_, * this->curvatureBilinearForm_);
//limitOutputFiles//      assembler.assemble (this->curvatureTgDiv_, * this->curvatureLinearForm_);
//limitOutputFiles//      assembler.assemble (* auxiliary_.vector(), * this->curvatureAdditionalForm_);
//limitOutputFiles//      dolfin::Vector curvatureRhs (this->curvatureTgDiv_);
//limitOutputFiles//      auxiliary_.vector()->get_local (addVals_, 2, addIdxs_);
//limitOutputFiles//      curvatureRhs.zero();
//limitOutputFiles//      curvatureRhs.set_local (addVals_, 2, addIdxs_);
//limitOutputFiles//      curvatureRhs += curvatureTgDiv_;
//limitOutputFiles//
//limitOutputFiles//      linearSolver_.solve (this->curvatureMatrix_, * auxiliary_.vector(), curvatureRhs);
//limitOutputFiles//
//limitOutputFiles//      print2csv (auxiliary_, problemData.savepath+"curvature.csv."+std::to_string(timeStep),
//limitOutputFiles//           pb_.meshManager_->orderedUXDofs().begin(),
//limitOutputFiles//           pb_.meshManager_->orderedUYDofs().begin(),
//limitOutputFiles//           pb_.meshManager_->orderedUXDofs().end(),
//limitOutputFiles//           auxiliary_.function_space()->dofmap()->tabulate_all_coordinates (* pb_.meshManager_->mesh())
//limitOutputFiles//           );
    }


    PostProcessor::~PostProcessor ()
    {
//addNotInMovingAbstract//      delete [] addVals_;
//addNotInMovingAbstract//      delete [] addIdxs_;
      delete [] uxDofsVals_; delete [] uyDofsVals_;
      delete [] pDofsVals_;
      delete [] wxDofsVals_; delete [] wyDofsVals_;
// The following deletes are NOT needed because they point to vectors managed by MeshManager
/*      delete [] uxDofsIdxs_; delete [] uyDofsIdxs_;
      delete [] pDofsIdxs_;
      delete [] wxDofsIdxs_; delete [] wyDofsIdxs_;
*/
    }

}//end of namespace
