#include "MovingTimeDependentProblem.h"
#include "discreteCurvature.h"
#include "discreteCurvature_onlyTP.h"
#include "myNavierstokesTimeCurvLinear.h"
#include "myNavierstokesTimeCurvLinearPreviousDomain.h"
#include "computeFreeSurfaceStress_onlyTP.h"
#include "utilities.h"

//#define GCL
  // if defined, assembles and file-writes the GCL dofs
  // NB: the initial files could be created anyway

namespace Ivan
{

    PostProcessor::PostProcessor (MovingTimeDependentProblem & pb) :
      PostProcessor::PostProcessor (pb, problemData.savepath+"mainEnergyForms.csv", problemData.savepath+"intGCL.csv", problemData.savepath+"divGCL.csv", problemData.savepath+"selectedDofs.csv")
    {}

    PostProcessor::PostProcessor (MovingTimeDependentProblem & pb, std::string outputFileName, std::string intGCLFileName, std::string divGCLFileName, std::string selectedDofsFileName) :
      pb_ (pb),
      coefficients_ (),
      formsNvalues_ (),
      formsNvaluesOnOld_ (),
      formsDepOnSol_ (),
      formsDepOnOld_ (),
      formsDepOnDispl_ (),
      formsDepOnDir_ (),
      gclSpace_ (* pb_.mesh()),
      gclForms_ (),
      outputFileName_ (outputFileName),
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
      additionalFormDofs_ (std::static_pointer_cast<Ivan::MovingLinearProblem 
                    < myNavierstokesTimeCurvLinear::FunctionSpace, computeFreeSurfaceStress_onlyTP::FunctionSpace,
                      myNavierstokesTimeCurvLinear::BilinearForm, myNavierstokesTimeCurvLinear::LinearForm,
                      myNavierstokesTimeCurvLinearPreviousDomain::LinearForm, computeFreeSurfaceStress_onlyTP::LinearForm > >
                    (pb_.timeSteppingProblem_)->additionalFormDofs ()),
      addIdxs_ (new dolfin::la_index [2]),
      addIdxs_w_ (new dolfin::la_index [2]),
      addVals_ (new double [2]),
      addVals_w_ (new double [2]),
      auxiliary_ (pb_.functionSpace())
    {
      // coefficients storing
      std::shared_ptr<dolfin::GenericFunction> coeffPtr;
      coeffPtr.reset (new dolfin::Constant (problemData.dt));
      coefficients_.emplace ("dt", coeffPtr);
      coeffPtr.reset (new dolfin::Constant (problemData.nu));
      coefficients_.emplace ("nu", coeffPtr);
      coeffPtr.reset (new dolfin::Constant (0,-inputData("gravityOn",1.0)));
      coefficients_.emplace ("gravityVector", coeffPtr);
      coeffPtr.reset (new dolfin::Constant (inputData("beta",0.0)));
      coefficients_.emplace ("beta", coeffPtr);
      coeffPtr.reset (new dolfin::Constant (0,inputData("wallVelY",0.0)));
      coefficients_.emplace ("wallVelocity", coeffPtr);
      coeffPtr.reset (new dolfin::Constant (inputData("gamma",0.0)));
      coefficients_.emplace ("gamma", coeffPtr);
      coeffPtr.reset (new dolfin::Constant (cos(inputData("thetaS",90.0)*3.14159265/180.0)*problemData.lx));
      //coeffPtr.reset (new dolfin::Constant (inputData("cosThetaS",cos(problemData.thetaS*3.14159265/180.0)*problemData.lx)));
      coefficients_.emplace ("cosThetaS", coeffPtr);
      coeffPtr.reset (new dolfin::Constant (inputData("stabPress",0.0)));
      coefficients_.emplace ("stabPress", coeffPtr);

      // forms construction and storing
      std::shared_ptr<dolfin::Form> enPtr;

      enPtr.reset (new mainEnergyForms::Form_kinEn (* pb_.mesh()));
      formsNvalues_.emplace ("kinEn", std::make_pair (enPtr, double ()));
      formsNvalues_["kinEn"].first->set_coefficient ("dt", dolfin::reference_to_no_delete_pointer (* coefficients_["dt"]));
      formsDepOnSol_.push_back (formsNvalues_["kinEn"].first);

      enPtr.reset (new mainEnergyForms::Form_grav (* pb_.mesh()));
      formsNvalues_.emplace ("grav", std::make_pair (enPtr, double ()));
      formsNvalues_["grav"].first->set_coefficient ("dt", dolfin::reference_to_no_delete_pointer (* coefficients_["dt"]));
      formsNvalues_["grav"].first->set_coefficient ("f", dolfin::reference_to_no_delete_pointer (* coefficients_["gravityVector"])); 

      enPtr.reset (new mainEnergyForms::Form_gammaLength (* pb_.mesh()));
      formsNvalues_.emplace ("gammaLength", std::make_pair (enPtr, double ()));
      formsNvalues_["gammaLength"].first->set_coefficient ("dt", dolfin::reference_to_no_delete_pointer (* coefficients_["dt"])); 
      formsNvalues_["gammaLength"].first->set_coefficient ("gamma", dolfin::reference_to_no_delete_pointer (* coefficients_["gamma"])); 

      enPtr.reset (new mainEnergyForms::Form_sigmaLength (* pb_.mesh()));
      formsNvalues_.emplace ("sigmaLength", std::make_pair (enPtr, double ()));
      formsNvalues_["sigmaLength"].first->set_coefficient ("dt", dolfin::reference_to_no_delete_pointer (* coefficients_["dt"])); 
      formsNvalues_["sigmaLength"].first->set_coefficient ("gamma", dolfin::reference_to_no_delete_pointer (* coefficients_["gamma"])); 
      formsNvalues_["sigmaLength"].first->set_coefficient ("cosThetaS", dolfin::reference_to_no_delete_pointer (* coefficients_["cosThetaS"])); 

      enPtr.reset (new mainEnergyForms::Form_viscPow (* pb_.mesh()));
      formsNvalues_.emplace ("viscPow", std::make_pair (enPtr, double ()));
      formsNvalues_["viscPow"].first->set_coefficient ("nu", dolfin::reference_to_no_delete_pointer (* coefficients_["nu"]));
      formsDepOnSol_.push_back (formsNvalues_["viscPow"].first);

      enPtr.reset (new mainEnergyForms::Form_nbc (* pb_.mesh()));
      formsNvalues_.emplace ("nbc", std::make_pair (enPtr, double ()));
      formsNvalues_["nbc"].first->set_coefficient ("beta", dolfin::reference_to_no_delete_pointer (* coefficients_["beta"])); 
      formsNvalues_["nbc"].first->set_coefficient ("wallVelocity", dolfin::reference_to_no_delete_pointer (* coefficients_["wallVelocity"])); 
      formsDepOnSol_.push_back (formsNvalues_["nbc"].first);

      enPtr.reset (new mainEnergyForms::Form_discr (* pb_.mesh()));
      formsNvalues_.emplace ("discr", std::make_pair (enPtr, double ()));
      formsNvalues_["discr"].first->set_coefficient ("dt", dolfin::reference_to_no_delete_pointer (* coefficients_["dt"]));
      formsDepOnSol_.push_back (formsNvalues_["discr"].first);
      formsDepOnOld_.push_back (formsNvalues_["discr"].first);

      enPtr.reset (new mainEnergyForms::Form_discrDiv (* pb_.mesh()));
      formsNvalues_.emplace ("discrDiv", std::make_pair (enPtr, double ()));
      formsNvalues_["discrDiv"].first->set_coefficient ("dt", dolfin::reference_to_no_delete_pointer (* coefficients_["dt"]));
      formsDepOnSol_.push_back (formsNvalues_["discrDiv"].first);
      formsDepOnOld_.push_back (formsNvalues_["discrDiv"].first);
      formsDepOnDispl_.push_back (formsNvalues_["discrDiv"].first);

      enPtr.reset (new mainEnergyForms::Form_epsg (* pb_.mesh()));
      formsNvaluesOnOld_.emplace ("epsg", std::make_pair (enPtr, double ()));
      formsNvaluesOnOld_["epsg"].first->set_coefficient ("dt", dolfin::reference_to_no_delete_pointer (* coefficients_["dt"]));
      formsNvaluesOnOld_["epsg"].first->set_coefficient ("f", dolfin::reference_to_no_delete_pointer (* coefficients_["gravityVector"])); 
      formsDepOnDispl_.push_back (formsNvaluesOnOld_["epsg"].first);

      /*enPtr.reset (new mainEnergyForms::Form_epsg (* pb_.mesh()));
      formsNvalues_.emplace ("epsg", std::make_pair (enPtr, double ()));
      formsNvalues_["epsg"].first->set_coefficient ("dt", dolfin::reference_to_no_delete_pointer (* coefficients_["dt"]));
      formsNvalues_["epsg"].first->set_coefficient ("f", dolfin::reference_to_no_delete_pointer (* coefficients_["gravityVector"])); 
      formsDepOnDispl_.push_back (formsNvalues_["epsg"].first);*/

      enPtr.reset (new mainEnergyForms::Form_tgDivw (* pb_.mesh()));
      formsNvaluesOnOld_.emplace ("tgDivw", std::make_pair (enPtr, double ()));
      formsNvaluesOnOld_["tgDivw"].first->set_coefficient ("gamma", dolfin::reference_to_no_delete_pointer (* coefficients_["gamma"])); 
      formsDepOnDispl_.push_back (formsNvaluesOnOld_["tgDivw"].first);

      enPtr.reset (new mainEnergyForms::Form_tgDivSol (* pb_.mesh()));
      formsNvalues_.emplace ("tgDivSol", std::make_pair (enPtr, double ()));
      formsNvalues_["tgDivSol"].first->set_coefficient ("gamma", dolfin::reference_to_no_delete_pointer (* coefficients_["gamma"])); 
      formsDepOnSol_.push_back (formsNvalues_["tgDivSol"].first);

      formsNvalues_.emplace ("supg", std::make_pair (enPtr, double ()));
      formsNvalues_["supg"].first->set_coefficient ("stabPress", dolfin::reference_to_no_delete_pointer (* coefficients_["stabPress"])); 
      formsDepOnSol_.push_back (formsNvalues_["supg"].first);

#ifdef GCL
      enPtr.reset (new mainEnergyForms::Form_intGCL (gclSpace_));
      gclForms_.emplace ("intGCL", enPtr);

      enPtr.reset (new mainEnergyForms::Form_divGCL (gclSpace_));
      gclForms_.emplace ("divGCL", enPtr);
      formsDepOnDispl_.push_back (gclForms_["divGCL"]);
#endif

      // Dirichlet datum and meshFunction setting
      std::shared_ptr<const dolfin::MeshFunction<std::size_t> > meshFunction (pb_.meshManager_->meshFacets());
      for (auto & formNvalue : formsNvalues_)
        formNvalue.second.first -> set_exterior_facet_domains (meshFunction);
      for (auto & formNvalue : formsNvaluesOnOld_)
        formNvalue.second.first -> set_exterior_facet_domains (meshFunction);

      std::ofstream balanceFile (balanceFileName_, std::ofstream::out);
      balanceFile << "kinEn, grav, gammaLength, sigmaLength, tpEnPot, kinEnFlux, gravFlux, viscPow, nbc, discr, discrDiv, epsg, tgDivw, tgDivSol, supg, epsGamma, epsPartGamma, Phi";
      balanceFile << std::endl;
      balanceFile.close();

      uxDofsIdxs_ = pb_.meshManager_->selectedOrderedUXDofs().data(); uyDofsIdxs_ = pb_.meshManager_->selectedOrderedUYDofs().data();
      pDofsIdxs_ = pb_.meshManager_->selectedOrderedPressDofs().data();
      wxDofsIdxs_ = pb_.meshManager_->selectedOrderedWXDofs().data(); wyDofsIdxs_ = pb_.meshManager_->selectedOrderedWYDofs().data();

const dolfin::GenericDofMap& dofmap_w (* pb_.meshManager_->displacement()->function_space()->dofmap());
dolfin::la_index numDofs_w (dofmap_w.dofs().size());
std::vector<dolfin::la_index> displacementDofs (dofmap_w.dofs());
std::vector<double> dofsCoords_w (dofmap_w.tabulate_all_coordinates (* pb_.meshManager_->mesh()));
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
      for (std::size_t i (0); i<2; ++i)
      {
        addIdxs_[i] = additionalFormDofs_[i];
        addIdxs_w_[i] = triplePointDofs_w[i];
      }
std::cerr << "PostProcessor::addIdxs_:addIdxs_w_ "; for (std::size_t iii(0); iii<2; ++iii) std::cerr << addIdxs_[iii] << ':' << addIdxs_w_[iii] << ' '; std::cerr << std::endl;
    }


    void PostProcessor::onOldDomain (int timeStep, const MovingTimeDependentProblem * const pb)
    {
      std::vector<double> bCoords (gclSpace_.dofmap()->tabulate_all_coordinates (* pb_.mesh()));
      std::vector<std::size_t> bIdxs (bCoords.size()/2);
      std::iota (bIdxs.begin(), bIdxs.end(), 0);
      dolfin::Vector bdiv;
#ifdef GCL
      gclForms_["divGCL"] -> set_coefficient ("w", dolfin::reference_to_no_delete_pointer (* pb_.meshManager_->displacement()));
      dolfin::assemble (bdiv, * gclForms_["divGCL"]);
      dolfin::Function bdivFun (gclSpace_);
      * bdivFun.vector() = bdiv;
      print2csv (bdivFun, divGCLFileName_+"old."+std::to_string(timeStep), bIdxs.begin(), bIdxs.end(), bCoords);
#endif

      const dolfin::Function & sol (pb_.solution()),
                          & oldSol (pb_.solution_.size()>1 ? pb_.solution_[pb_.solution_.size()-2].second : sol),
                               & w (* pb_.meshManager_->displacement());
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
    }

    void PostProcessor::operator() (int timeStep, const MovingTimeDependentProblem * const pb)
    {
      const dolfin::Function & sol (pb_.solution()),
                          & oldSol (pb_.solution_.size()>1 ? pb_.solution_[pb_.solution_.size()-2].second : sol),
                               & w (* pb_.meshManager_->displacement());
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
//        std::cerr << "processing " << formNvalue.first << std::endl;
        formNvalue.second.second = dolfin::assemble (* formNvalue.second.first);
      }

      std::ofstream outputFile (outputFileName_, std::ofstream::app);
      balanceFile.setf ( std::ios::scientific, std::ios::floatfield);
      balanceFile.precision (15);
//      balanceFile << "kinEn, grav, gammaLength, sigmaLength, tpEnPot, kinEnFlux, gravFlux, viscPow, nbc, discr, discrDiv, epsg, tgDivw, tgDivSol, supg, epsGamma, epsPartGamma, Phi";
      balanceFile << formsNvalues_["kinEn"].second << ',';
      balanceFile << formsNvalues_["grav"].second << ',';
      balanceFile << formsNvalues_["gammaLength"].second << ',';
      balanceFile << formsNvalues_["sigmaLength"].second << ',';
      balanceFile << problemData.gamma*cos(problemData.thetaS*3.14159265/180.0)*problemData.lx*w;
      balanceFile << formsNvalues_["kinEnFlux"].second << ',';
      balanceFile << formsNvalues_["gravFlux"].second << ',';
      balanceFile << formsNvalues_["viscPow"].second << ',';
      balanceFile << formsNvalues_["nbc"].second << ',';
      outputFile << formsNvalues_["discr"].second << ',';
      outputFile << formsNvalues_["discrDiv"].second << ',';
      outputFile << formsNvalues_["bdTerm"].second << ',';
      outputFile << formsNvalues_["enPot"].second << ',';
      outputFile << formsNvaluesOnOld_["epsg"].second << ',';
      outputFile << formsNvalues_["epsg"].second << ',';
      outputFile << formsNvalues_["bdGravu"].second <<  ',';
      outputFile << formsNvaluesOnOld_["bdGravw"].second <<  ',';
      outputFile << formsNvalues_["bdGravw"].second <<  ',';
      outputFile << formsNvaluesOnOld_["divEpsGamma"].second << ',';
      outputFile << (formsNvaluesOnOld_["divEpsGamma"].second = dolfin::assemble (* formsNvaluesOnOld_["divEpsGamma"].first)) << ',';
      outputFile << formsNvaluesOnOld_["divEpsGammaErr"].second << ',';
      outputFile << formsNvalues_["tgDivSol"].second << ',';
      outputFile << formsNvalues_["inflow"].second << ',';
      outputFile << formsNvalues_["supg"].second << std::endl;
      outputFile.close();

      std::vector<double> bCoords (gclSpace_.dofmap()->tabulate_all_coordinates (* pb_.mesh()));
      std::vector<std::size_t> bIdxs (bCoords.size()/2);
      std::iota (bIdxs.begin(), bIdxs.end(), 0);
#ifdef GCL
      dolfin::Vector b;
      dolfin::assemble (b, * gclForms_["intGCL"]);
      dolfin::Function bFun (gclSpace_);
      * bFun.vector() = b;
      print2csv (bFun, intGCLFileName_+"."+std::to_string(timeStep), bIdxs.begin(), bIdxs.end(), bCoords);
      dolfin::Vector bdiv;
      dolfin::assemble (bdiv, * gclForms_["divGCL"]);
      * bFun.vector() = bdiv;
      print2csv (bFun, divGCLFileName_+"."+std::to_string(timeStep), bIdxs.begin(), bIdxs.end(), bCoords);
#endif

std::cerr<<"uxDofsIdxs ("<<uxDofsNum_<< ") "; for (dolfin::la_index i (0); i<uxDofsNum_; ++i) std::cerr << uxDofsIdxs_[i] << ", "; std::cerr<<std::endl;
std::cerr<<"uyDofsIdxs ("<<uyDofsNum_<< ") "; for (dolfin::la_index i (0); i<uyDofsNum_; ++i) std::cerr << uyDofsIdxs_[i] << ", "; std::cerr<<std::endl;
std::cerr<<"pDofsIdxs ("<<pDofsNum_<< ") "; for (dolfin::la_index i (0); i<pDofsNum_; ++i) std::cerr << pDofsIdxs_[i] << ", "; std::cerr<<std::endl;
std::cerr<<"wxDofsIdxs ("<<wxDofsNum_<< ") "; for (dolfin::la_index i (0); i<wxDofsNum_; ++i) std::cerr << wxDofsIdxs_[i] << ", "; std::cerr<<std::endl;
std::cerr<<"wyDofsIdxs ("<<wyDofsNum_<< ") "; for (dolfin::la_index i (0); i<wyDofsNum_; ++i) std::cerr << wyDofsIdxs_[i] << ", "; std::cerr<<std::endl;
      sol.vector()->get (uxDofsVals_, & uxDofsNum_, & uxDofsIdxs_);
      sol.vector()->get (uyDofsVals_, & uyDofsNum_, & uyDofsIdxs_);
      sol.vector()->get (pDofsVals_, & pDofsNum_, & pDofsIdxs_);
      w.vector()->get (wxDofsVals_, & wxDofsNum_, & wxDofsIdxs_);
      w.vector()->get (wyDofsVals_, & wyDofsNum_, & wyDofsIdxs_);
      print2csv (std::vector<const double *> ({uxDofsVals_, uyDofsVals_, pDofsVals_, wxDofsVals_, wyDofsVals_}),
                 std::vector<const dolfin::la_index *> ({& uxDofsNum_, & uyDofsNum_, & pDofsNum_, & wxDofsNum_, & wyDofsNum_}),
                 selectedDofsFileName_+"."+std::to_string(timeStep));

      //discrete curvature

      dolfin::Assembler assembler;
std::cerr << " OK fino a " << __LINE__ << std::endl;
      assembler.assemble (this->curvatureMatrix_, * this->curvatureBilinearForm_);
std::cerr << " OK fino a " << __LINE__ << std::endl;
      assembler.assemble (this->curvatureTgDiv_, * this->curvatureLinearForm_);
std::cerr << " OK fino a " << __LINE__ << std::endl;
      assembler.assemble (* auxiliary_.vector(), * this->curvatureAdditionalForm_);
std::cerr << " OK fino a " << __LINE__ << std::endl;
      dolfin::Vector curvatureRhs (this->curvatureTgDiv_);
      auxiliary_.vector()->get_local (addVals_, 2, addIdxs_);
      curvatureRhs.zero();
      curvatureRhs.set_local (addVals_, 2, addIdxs_);
      curvatureRhs += curvatureTgDiv_;

      linearSolver_.solve (this->curvatureMatrix_, * auxiliary_.vector(), curvatureRhs);

      print2csv (auxiliary_, problemData.savepath+"curvature.csv."+std::to_string(timeStep),
           pb_.meshManager_->orderedUXDofs().begin(),
           pb_.meshManager_->orderedUYDofs().begin(),
           pb_.meshManager_->orderedUXDofs().end(),
           auxiliary_.function_space()->dofmap()->tabulate_all_coordinates (* pb_.meshManager_->mesh())
           );
    }


    PostProcessor::~PostProcessor ()
    {
      delete [] addVals_;
      delete [] addIdxs_;
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
