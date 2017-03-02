#ifndef MESHMANAGER_H_INCLUDE_GUARD
#define MESHMANAGER_H_INCLUDE_GUARD

#include <dcp/differential_problems/AbstractProblem.h>
#include <dcp/differential_problems/LinearProblem.h>
#include <dcp/differential_problems/utilities.h> //#include "utilities.h"
#include <dcp/differential_problems/meshManagerUflTools.h>
#include <limits>
#include <algorithm>
#include <iterator>
#include <string>

template <class T_MeshMover> class DepthEvaluator; //fwd declaration
extern struct ProblemData problemData;

namespace dcp
{

template <class T_MeshMover = dolfin::ALE>
class MeshManager
{

  public:

    typedef T_MeshMover MeshMover;

    MeshManager () = delete;
  
  	//! Constructor
  	/*! 
  	 *  \param mesh the mesh to be handled
  	 *	\param meshFacets the MeshFunction labelling the mesh, with the following values:
 	 *		1 : the boundary on which the given displacement has to be imposed
	 * 		2 : fixed boundary
	 *		3 : slip boundary
	 *		0 : interior facets / never mind
	 *	\param noSlipComponent the function component along which slip is not allowed,
	 *	on the boundary marked as 3 by meshFacets ( default value = std::numeric_limits<std::size_t>::max() )
	 */
  	MeshManager (const std::shared_ptr<dolfin::Mesh> mesh,
				 const std::shared_ptr<const dolfin::FacetFunction<std::size_t> > meshFacets,
         const std::shared_ptr<dcp::AbstractProblem> aleProblem,
         const dolfin::FunctionSpace displCoeffSpace,
				 const std::size_t noSlipComponent=std::numeric_limits<std::size_t>::max());
    //?? noSlipComponet

    //! DEPRECATED (use #storeOrderedDofIdxs, instead)
    void getDofOrder (const dolfin::FunctionSpace & funSp, std::size_t component=0);
    //! Initializer ...??
    void storeOrderedDofIdxs (const dolfin::FunctionSpace & funSp, std::list<std::string> names, int component=-1);
    void recursiveStoreOrderedDofIdxs (const dolfin::FunctionSpace & funSp, std::list<std::string> & names, const std::vector<double> & dofsCoords, std::map<std::string, std::vector<dolfin::la_index> > & orderedDofs);

    /*************************** GETTERS ***********************/

    //!@name getters of DEPRECATED members
    //!@{
    const std::vector<std::size_t> & orderedDisplacementDofs () const
    {
        return orderedDisplacementDofs_;
    }
    const std::vector<std::pair<double, double> > & orderedDisplacementDofsCoords () const
    {
        return orderedDisplacementDofsCoords_;
    }
    const std::vector<std::size_t> & orderedUXDofs () const { return orderedUXDofs_; }
    const std::vector<std::size_t> & orderedUYDofs () const { return orderedUYDofs_; }
    const std::vector<std::size_t> & orderedPressDofs () const { return orderedPressDofs_; }
	  const std::vector<std::size_t> & orderedWXDofs () const { return orderedWXDofs_; }
    const std::vector<std::size_t> & orderedWYDofs () const { return orderedWYDofs_; }
    const std::vector<dolfin::la_index> & selectedOrderedUXDofs () const { return selectedOrderedUXDofs_; }
    const std::vector<dolfin::la_index> & selectedOrderedUYDofs () const { return selectedOrderedUYDofs_; }
    const std::vector<dolfin::la_index> & selectedOrderedPressDofs () const { return selectedOrderedPressDofs_; }
	  const std::vector<dolfin::la_index> & selectedOrderedWXDofs () const { return selectedOrderedWXDofs_; }
    const std::vector<dolfin::la_index> & selectedOrderedWYDofs () const { return selectedOrderedWYDofs_; }
    //!@}
    
    const std::map<std::string, std::vector<dolfin::la_index> > & orderedDofs () const { return orderedDofs_; }
    const std::map<std::string, std::vector<dolfin::la_index> > & orderedDisplDofs () const { return orderedDisplDofs_; }

//    const dolfin::FunctionSpace& functionSpace() const;
    std::shared_ptr<const dolfin::Mesh> mesh() const;
    std::shared_ptr<const dolfin::Mesh> oldMesh() const;
    std::shared_ptr<const dolfin::MeshFunction<std::size_t> > meshFacets() const;
	
//    void moveMesh (const std::string& component="all")
//    {
//        moveMesh (this->w_, component);
//    }
//    void moveMesh (const dolfin::Function& displacement, const std::string& component="all", const double dt=1.0);
    std::shared_ptr <const dolfin::Function> moveMesh (const std::string& component="all")
    {
        return moveMesh (this->w_, component);
    }
    std::shared_ptr <const dolfin::Function> moveMesh (const dolfin::Function& displacement, const std::string& component="all", const double dt=1.0);
    std::shared_ptr <const dolfin::Function> computeDisplacement (const dolfin::Function& displacement, const std::string& component="all", const double dt=1.0);

    std::shared_ptr <const dolfin::Function> displacement () const
	{ return dolfin::reference_to_no_delete_pointer (w_); }
    
    void setImposedDisplacement (std::shared_ptr<const dolfin::FunctionSpace> imposedDisplacementSpace);

    //! Preserves previous labelling
    void initializeMesh (const dolfin::Mesh & mesh);

    //! Get displacement function space
    const dolfin::FunctionSpace & displacementFunctionSpace () const;

    //! Set displacement
    /*!
     *  Assumes the displacement given in input lives in the FunctionSpace given by this->displacementFunctionSpace()
     *  \param displacement
     */
    void setDisplacement (const dolfin::Function & displacement);

  	std::shared_ptr<const dolfin::FacetFunction<std::size_t> > meshFacets () { return meshFacets_; }

    void print2csv (const dolfin::Function & fun, const std::string & prependToFileName, const std::string & appendToFileName, bool printAlsoDisplacement) const;
    void print2csv (const dolfin::Function & fun, const std::string & prependToFileName, const std::string & appendToFileName, const std::map<std::string,std::vector<dolfin::la_index> > & dofs) const;

  protected:

    T_MeshMover meshMover_;
    const std::shared_ptr<dolfin::Mesh> mesh_;
    meshManagerUflTools::FunctionSpace wFunSp_; //dolfin::FunctionSpace wFunSp_;
    dolfin::Function w_;
    const dolfin::FunctionSpace & displCoeffSpace_;
    const std::shared_ptr <dcp::AbstractProblem> aleProblem_;
//REM//    dcp::LinearProblem<laplaceALE::BilinearForm, laplaceALE::LinearForm> laplacePb_;
    bool displacementIsComputed_;
//    std::shared_ptr<dcp::AbstractProblem> problemALE_;
/*    std::vector<std::shared_ptr<dolfin::FiniteElement> > subElements_;
    std::vector<std::shared_ptr<dolfin::GenericDofMap> > subDofmaps_;
*/
  	const std::shared_ptr<const dolfin::FacetFunction<std::size_t> > meshFacets_;
  	const std::size_t noSlipComponent_;
    dolfin::MeshFunction<bool> on_boundary_;
    dolfin::MeshFunction<bool> on_to_be_neglected_boundary_;

    //!@name DEPRECATED members
    //!@{
    std::vector<std::size_t> orderedMeshIdxs_;
    std::vector<std::size_t> orderedDisplacementDofs_;
    std::vector<std::pair<double, double> > orderedDisplacementDofsCoords_;
    std::vector<std::size_t> orderedUXDofs_;
    std::vector<std::size_t> orderedUYDofs_;
    std::vector<std::size_t> orderedPressDofs_;
    std::vector<std::size_t> orderedWXDofs_;
    std::vector<std::size_t> orderedWYDofs_;
    std::vector<dolfin::la_index> selectedOrderedUXDofs_;
    std::vector<dolfin::la_index> selectedOrderedUYDofs_;
    std::vector<dolfin::la_index> selectedOrderedPressDofs_;
    std::vector<dolfin::la_index> selectedOrderedWXDofs_;
    std::vector<dolfin::la_index> selectedOrderedWYDofs_;
    //!@}

 std::map<std::string,std::vector<dolfin::la_index> > orderedDofs_;
 std::map<std::string,std::vector<dolfin::la_index> > orderedDisplDofs_;

    //ImposedDisplacement imposedDisplacement_;
    //ImposedFromFile imposedDisplacement_;
    ImposedDisplacement imposedDisplacement_;
    std::shared_ptr<const dolfin::FunctionSpace> imposedDisplacementSpace_;
    bool displacementIsImposed_;

    friend class DepthEvaluator<T_MeshMover>;

    dolfin::Mesh & oldMesh_;

    dolfin::Constant uno_, zero_, penalty_;
};

template <class T_MeshMover = dolfin::ALE>
class DepthEvaluator
{

    public:
        DepthEvaluator () = delete;
        DepthEvaluator (const MeshManager<T_MeshMover> * mm) :
          meshManager_ (mm)
          {}

        void operator() (dolfin::Array<double>& values, const dolfin::Array<double>& x, const double& t)
        {
            values[0] = 0;
            values[1] = - 9.81 * (std::find_if
                                  ( meshManager_->orderedDisplacementDofsCoords_.cbegin(), meshManager_->orderedDisplacementDofsCoords_.cend(),
                                    [&x, &problemData] (std::pair<double,double> xy)
                                    {
                                      return dolfin::near (x[0], xy.first, 0.1*problemData.lx/problemData.nx);
                                    }
                                  ) -> second) ;
        }

    protected:
        const MeshManager<T_MeshMover> * meshManager_;
};/*class DepthEvaluator : public dolfin::Expression
{
public:

  void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x) const
  {
      values[0] = 0;
      values[1] = - 9.81 * find_if
                           ( meshManager_->orderedDisplacementDofsCoords_.cbegin(), meshManager_->orderedDisplacementDofsCoords_.cend(),
                             [&x, &problemData] (std::pair<double,double> xy)
                             {
                               return dolfin::near (x[0], xy.first, 0.1*problemData.lx/problemData.nx);
                             }
                           ) -> second;
  }

  std::size_t value_rank() const
  {
    return 1;
  }

  std::size_t value_dimension(std::size_t i) const
  {
    return 2;
  }

  void setMeshManager (MeshManager<T_MeshMover> & mm)
  {
    this->meshManager_.reset (mm);
  }

protected:
  
  std::shared_ptr<MeshManager<T_MeshMover> > meshManager_;
  
};*/

// =============== IMPLEMENTATION ===============

template <class T_MeshMover>
MeshManager<T_MeshMover>::MeshManager (const std::shared_ptr<dolfin::Mesh> mesh,
 									   const std::shared_ptr<const dolfin::FacetFunction<std::size_t> > meshFacets,
                     const std::shared_ptr<dcp::AbstractProblem> aleProblem,
                     const dolfin::FunctionSpace displCoeffSpace,
									   const std::size_t noSlipComponent) :
  meshMover_ (),
	mesh_ (mesh),
	//wFunSp_ (laplaceVec::FunctionSpace(*mesh_)),
	wFunSp_ (mesh_), //SPACE//displCoeffSpace), //myNavierstokesTimeCurvLinear::CoefficientSpace_w (*mesh_)),
	w_ (wFunSp_),
  displCoeffSpace_ (displCoeffSpace),
//REM//	aleFunSp_ (laplaceALE::FunctionSpace (*mesh_)),
  aleProblem_ (aleProblem),
//REM//  laplacePb_ (laplaceFunSp_),
  displacementIsComputed_ (false),
//	  problemALE_ (new dcp::LinearProblem<laplaceVec::BilinearForm, laplaceVec::LinearForm> (dolfin::reference_to_no_delete_pointer(wFunSp_))),
	meshFacets_ (meshFacets),
	noSlipComponent_ (noSlipComponent),
	on_boundary_ (* mesh_, 0),
	on_to_be_neglected_boundary_ (* mesh_, 0),
    orderedMeshIdxs_ (mesh_->coordinates().size()/2),
    orderedDisplacementDofs_ (),
    orderedDisplacementDofsCoords_ (),
    orderedUXDofs_ (),
    orderedUYDofs_ (),
    orderedPressDofs_ (),
    orderedWXDofs_ (),
    orderedWYDofs_ (),
    selectedOrderedUXDofs_ (),
    selectedOrderedUYDofs_ (),
    selectedOrderedPressDofs_ (),
    selectedOrderedWXDofs_ (),
    selectedOrderedWYDofs_ (),
 orderedDofs_(),
 orderedDisplDofs_(),
    displacementIsImposed_ (false),
    oldMesh_ (* mesh),
    uno_ (1.0),
    zero_ (0.0),
    penalty_ (1.0e20)
{
std::cerr << displCoeffSpace_.str(true) << std::endl;
    const std::vector<double> & coords (mesh_->coordinates());
//for (auto it=coords.begin(); it!=coords.end(); ++it) std::cerr << *it << ' '; std::cerr << std::endl;
    std::vector<std::tuple<std::size_t,double,double> > ordered;
    for (std::size_t k=0; k!=coords.size()/2; ++k)
      ordered.push_back (std::make_tuple (k, coords[2*k], coords[2*k+1]));
    std::sort (ordered.begin(),ordered.end(),compare);
    std::transform (ordered.begin(),ordered.end(),orderedMeshIdxs_.begin(),getFirst);

    const std::vector<double> displCoords (wFunSp_.dofmap()->tabulate_all_coordinates(* mesh_));
    std::list<std::string> names ({"wx","wy"});
    recursiveStoreOrderedDofIdxs (wFunSp_, names, displCoords, orderedDisplDofs_);

    dolfin::Constant zeroVec (0.0,0.0);
/*REM//    laplacePb_->setCoefficient ("bilinear_form",dolfin::reference_to_no_delete_pointer(uno_),"k");
    //laplacePb_->setCoefficient ("bilinear_form",dolfin::reference_to_no_delete_pointer(zero_),"k");
    laplacePb_->setCoefficient ("bilinear_form",dolfin::reference_to_no_delete_pointer(penalty_),"penalty");
    laplacePb_->setCoefficient ("linear_form",dolfin::reference_to_no_delete_pointer(uno_),"k");
    //laplacePb_->setCoefficient ("linear_form",dolfin::reference_to_no_delete_pointer(zero_),"k");
    laplacePb_->setCoefficient ("linear_form",dolfin::reference_to_no_delete_pointer(zero_),"f");
    laplacePb_->setCoefficient ("linear_form",dolfin::reference_to_no_delete_pointer(zero_),"g");
    laplacePb_->setCoefficient ("linear_form",dolfin::reference_to_no_delete_pointer(penalty_),"penalty");
    
		laplacePb_->setIntegrationSubdomain ("bilinear_form", meshFacets_, dcp::SubdomainType::BOUNDARY_FACETS);
		laplacePb_->setIntegrationSubdomain ("linear_form", meshFacets_, dcp::SubdomainType::BOUNDARY_FACETS);
	  BottomBd bottomBoundary;
    dolfin::DirichletBC dirBCbottom (laplaceFunSp_, zero_, bottomBoundary, "geometric");
std::unordered_map<std::size_t, double> bdval;
dirBCbottom.get_boundary_values(bdval);
    laplacePb_->addDirichletBC (dirBCbottom, std::string("bottom"));
*/
}

template <class T_MeshMover>
void MeshManager<T_MeshMover>::getDofOrder (const dolfin::FunctionSpace & funSp, std::size_t component)
{
    const std::vector<double> & dofsCoords (funSp.dofmap()->tabulate_all_coordinates (* mesh_));
    std::vector<dolfin::la_index> dofs (funSp[component]->dofmap()->dofs());
    std::vector<std::tuple<std::size_t,double,double,std::size_t> > ordered;
    const std::vector<double> & displDofsCoords (wFunSp_.dofmap()->tabulate_all_coordinates (* mesh_));
//std::cerr << "==========" << std::endl;
    for (std::size_t subcomp=0; subcomp!=funSp[component]->element()->value_dimension(0); ++subcomp)
    {
      std::vector<dolfin::la_index> dofs ((* funSp[component])[subcomp]->dofmap()->dofs());
      for (auto it=dofs.begin(); it!=dofs.end(); ++it)
      {
        ordered.push_back (std::make_tuple (*it, dofsCoords[2*(*it)], dofsCoords[2*(*it)+1], subcomp));
//std::cerr << *it << ' ' << dofsCoords[2*(*it)] << ' ' << dofsCoords[2*(*it)+1] << ' ' << subcomp << std::endl;
      }
    }
    std::sort (ordered.begin(),ordered.end(),compareDofs);
    orderedDisplacementDofs_.resize (ordered.size());
    orderedDisplacementDofsCoords_.resize (ordered.size());
    std::transform (ordered.begin(),ordered.end(),orderedDisplacementDofs_.begin(),getFirstDofs);
    std::transform (ordered.begin(),ordered.end(),orderedDisplacementDofsCoords_.begin(),getDofsCoords);

    std::vector<std::tuple<std::size_t,double,double,std::size_t> > orderedUX;
    std::vector<dolfin::la_index> dofsUX ((* funSp[0])[0]->dofmap()->dofs());
    for (auto it (dofsUX.begin()); it!=dofsUX.end(); ++it)
        orderedUX.push_back (std::make_tuple (*it, dofsCoords[2*(*it)], dofsCoords[2*(*it)+1], 2));
//std::cerr << "==========" << std::endl;
    std::sort (orderedUX.begin(),orderedUX.end(),compareDofs);
    orderedUXDofs_.resize (orderedUX.size());
    std::transform (orderedUX.begin(),orderedUX.end(),orderedUXDofs_.begin(),getFirstDofs);

    std::vector<std::tuple<std::size_t,double,double,std::size_t> > orderedUY;
    std::vector<dolfin::la_index> dofsUY ((* funSp[0])[1]->dofmap()->dofs());
    for (auto it (dofsUY.begin()); it!=dofsUY.end(); ++it)
        orderedUY.push_back (std::make_tuple (*it, dofsCoords[2*(*it)], dofsCoords[2*(*it)+1], 2));
//std::cerr << "==========" << std::endl;
    std::sort (orderedUY.begin(),orderedUY.end(),compareDofs);
    orderedUYDofs_.resize (orderedUY.size());
    std::transform (orderedUY.begin(),orderedUY.end(),orderedUYDofs_.begin(),getFirstDofs);

    std::vector<std::tuple<std::size_t,double,double,std::size_t> > orderedPress;
    std::vector<dolfin::la_index> dofsPress (funSp[1]->dofmap()->dofs());
    for (auto it (dofsPress.begin()); it!=dofsPress.end(); ++it)
        orderedPress.push_back (std::make_tuple (*it, dofsCoords[2*(*it)], dofsCoords[2*(*it)+1], 2));
//std::cerr << "==========" << std::endl;
    std::sort (orderedPress.begin(),orderedPress.end(),compareDofs);
    orderedPressDofs_.resize (orderedPress.size());
    std::transform (orderedPress.begin(),orderedPress.end(),orderedPressDofs_.begin(),getFirstDofs);

    std::vector<std::tuple<std::size_t,double,double,std::size_t> > orderedWX;
    std::vector<dolfin::la_index> dofsWX (wFunSp_[0]->dofmap()->dofs());
    for (auto it (dofsWX.begin()); it!=dofsWX.end(); ++it)
        orderedWX.push_back (std::make_tuple (*it, displDofsCoords[2*(*it)], displDofsCoords[2*(*it)+1], 2));
//std::cerr << "==========" << std::endl;
    std::sort (orderedWX.begin(),orderedWX.end(),compareDofs);
    orderedWXDofs_.resize (orderedWX.size());
    std::transform (orderedWX.begin(),orderedWX.end(),orderedWXDofs_.begin(),getFirstDofs);

    std::vector<std::tuple<std::size_t,double,double,std::size_t> > orderedWY;
    std::vector<dolfin::la_index> dofsWY (wFunSp_[1]->dofmap()->dofs());
    for (auto it (dofsWY.begin()); it!=dofsWY.end(); ++it)
        orderedWY.push_back (std::make_tuple (*it, displDofsCoords[2*(*it)], displDofsCoords[2*(*it)+1], 2));
//std::cerr << "==========" << std::endl;
    std::sort (orderedWY.begin(),orderedWY.end(),compareDofs);
    orderedWYDofs_.resize (orderedWY.size());
    std::transform (orderedWY.begin(),orderedWY.end(),orderedWYDofs_.begin(),getFirstDofs);

//std::cerr << "orderedMeshIdxs_ (" << orderedMeshIdxs_.size() << ")   "; for (auto it=orderedMeshIdxs_.begin(); it!=orderedMeshIdxs_.end(); ++it) std::cerr << *it << ' '; std::cerr << std::endl;
//std::cerr << "orderedDisplacementDofs_ (" << orderedDisplacementDofs_.size() << ")   "; for (auto it=orderedDisplacementDofs_.begin(); it!=orderedDisplacementDofs_.end(); ++it) std::cerr << *it << ' '; std::cerr << std::endl;
//std::cerr << "dofsCoords (" << dofsCoords.size() << ")   "; for (auto it=dofsCoords.begin(); it!=dofsCoords.end(); ++it) std::cerr << *it << ' '; std::cerr << std::endl;
//std::cerr << " OK fino a " << __LINE__ << std::endl;
    std::ofstream dofFile (problemData.savepath+"dof_ordered.dat", std::ofstream::out);
    dofFile << dofsCoords.size() << ',';
    for (auto& d : dofsCoords)
      dofFile << d << ',';
    dofFile << std::endl << dofs.size() << ',';
    for (auto& d : dofs)
      dofFile << d << ',';
    dofFile << std::endl << orderedDisplacementDofs_.size() << ',';
    for (auto& d : orderedDisplacementDofs_)
      dofFile << d << ',';
    dofFile << std::endl << orderedDisplacementDofsCoords_.size() << ',';
    for (auto& d : orderedDisplacementDofsCoords_)
      dofFile << d.first << ',' << d.second << ",";
    dofFile << std::endl << orderedUXDofs_.size() << ',';
    for (auto& d : orderedUXDofs_)
      dofFile << d << ',';
    dofFile << std::endl << orderedUYDofs_.size() << ',';
    for (auto& d : orderedUYDofs_)
      dofFile << d << ',';
    dofFile << std::endl << orderedPressDofs_.size() << ',';
    for (auto& d : orderedPressDofs_)
      dofFile << d << ',';
    dofFile << std::endl << orderedWXDofs_.size() << ',';
    for (auto& d : orderedWXDofs_)
      dofFile << d << ',';
    dofFile << std::endl << orderedWYDofs_.size() << ',';
    for (auto& d : orderedWYDofs_)
      dofFile << d << ',';
    dofFile << std::endl;
    dofFile.close();

//!!! Seleziona solo i dof delle patch dei punti topleft, topmid, topright
    /*std::size_t topleft  (2*problemData.ny*(2*problemData.nx+1)),         topleftP1  (problemData.ny*(problemData.nx+1)),
                topright ((2*problemData.ny+1)*(2*problemData.nx+1)-1),   toprightP1 ((problemData.ny+1)*(problemData.nx+1)-1),
                topmid ((topleft+topright)/2),                          topmidP1   ((topleftP1+toprightP1)/2);*/
    // P1!
    std::size_t topleft  (problemData.ny*(problemData.nx+1)),         topleftP1  (problemData.ny*(problemData.nx+1)),
                topright ((problemData.ny+1)*(problemData.nx+1)-1),   toprightP1 ((problemData.ny+1)*(problemData.nx+1)-1),
                topmid ((topleft+topright)/2),                          topmidP1   ((topleftP1+toprightP1)/2);
    std::vector<std::size_t> tmp ({topleft,topleft+1,topleft+2, topmid-2,topmid-1,topmid,topmid+1,topmid+2, topright-2,topright-1,topright}),
                             tmpP1 ({topleftP1,topleftP1+1,     topmidP1-1,topmidP1,topmidP1+1,             toprightP1-1,toprightP1});
    std::vector<std::size_t> selectedDofPos, selectedDofPosP1;
    selectedDofPos.insert (selectedDofPos.end(), tmp.begin(), tmp.end());
    selectedDofPosP1.insert (selectedDofPosP1.end(), tmpP1.begin(), tmpP1.end());
    std::transform (tmp.begin(), tmp.end(), tmp.begin(), [&problemData] (std::size_t el) -> std::size_t { return el-(problemData.nx+1); } );
    std::transform (tmpP1.begin(), tmpP1.end(), tmpP1.begin(), [&problemData] (std::size_t el) { return el-problemData.nx-1; } );
    selectedDofPos.insert (selectedDofPos.end(), tmp.begin(), tmp.end());
    selectedDofPosP1.insert (selectedDofPosP1.end(), tmpP1.begin(), tmpP1.end());
    std::transform (tmp.begin(), tmp.end(), tmp.begin(), [&problemData] (std::size_t el) { return el-(problemData.nx+1); } );
    std::transform (tmpP1.begin(), tmpP1.end(), tmpP1.begin(), [&problemData] (std::size_t el) { return el-problemData.nx-1; } );
    selectedDofPos.insert (selectedDofPos.end(), tmp.begin(), tmp.end());
    selectedDofPosP1.insert (selectedDofPosP1.end(), tmpP1.begin(), tmpP1.end());
std::cerr << "selectedDofPos (" << selectedDofPos.size() << ") "; for (auto e : selectedDofPos) std::cerr << e << ", "; std::cerr << std::endl;
    
    for (std::size_t idx : selectedDofPos)
    {
      selectedOrderedUXDofs_.push_back (orderedUXDofs_[idx]);
      selectedOrderedUYDofs_.push_back (orderedUYDofs_[idx]);
    }
    for (std::size_t idx : selectedDofPosP1)
    {
      selectedOrderedPressDofs_.push_back (orderedPressDofs_[idx]);
      selectedOrderedWXDofs_.push_back (orderedWXDofs_[idx]);
      selectedOrderedWYDofs_.push_back (orderedWYDofs_[idx]);
    }
}

template <class T_MeshMover>
void MeshManager<T_MeshMover>::storeOrderedDofIdxs (const dolfin::FunctionSpace & funSp, std::list<std::string> names, int component)
{
    const std::vector<double> & dofsCoords (funSp.dofmap()->tabulate_all_coordinates (* mesh_));

    if (-1 == component)
      recursiveStoreOrderedDofIdxs (funSp, names, dofsCoords, orderedDofs_);
    else
      recursiveStoreOrderedDofIdxs (* funSp[component], names, dofsCoords, orderedDofs_);

    return;
}

template <class T_MeshMover>
void MeshManager<T_MeshMover>::recursiveStoreOrderedDofIdxs (const dolfin::FunctionSpace & funSp, std::list<std::string> & names, const std::vector<double> & dofsCoords, std::map<std::string, std::vector<dolfin::la_index> > & stored)
{
  std::size_t numComp (funSp.element()->num_sub_elements());

  if (0 == numComp)
  {
    std::vector<dolfin::la_index> dofs (funSp.dofmap()->dofs());
    std::vector<std::tuple<dolfin::la_index,double,double> > orderingTool(dofs.size());
    auto itT (orderingTool.begin());
    std::transform (dofs.begin(), dofs.end(), orderingTool.begin(), [& dofsCoords](dolfin::la_index idx)
      { return std::make_tuple (idx, dofsCoords[2*idx], dofsCoords[2*idx+1]); });
    std::sort (orderingTool.begin(),orderingTool.end(),compare);

    stored.insert (std::make_pair(names.front(), std::vector<dolfin::la_index> (dofs.size())));
    std::vector<dolfin::la_index> & ordDofs (stored[names.front()]);
    std::transform (orderingTool.begin(), orderingTool.end(), ordDofs.begin(), [](decltype(* orderingTool.begin()) & t)
      { return std::get<0> (t); } );

    names.pop_front();
    return;
  }

  for (std::size_t comp (0); comp!=numComp; ++comp)
    recursiveStoreOrderedDofIdxs (* funSp[comp], names, dofsCoords, stored);
}

template <class T_MeshMover>
std::shared_ptr<const dolfin::Mesh> MeshManager<T_MeshMover>::mesh() const
{
    return mesh_;
}
template <class T_MeshMover>
std::shared_ptr<const dolfin::Mesh> MeshManager<T_MeshMover>::oldMesh() const
{
    return dolfin::reference_to_no_delete_pointer (oldMesh_);
}
template <class T_MeshMover>
std::shared_ptr<const dolfin::MeshFunction<std::size_t> > MeshManager<T_MeshMover>::meshFacets() const
{
    return dolfin::reference_to_no_delete_pointer (meshFacets_);
}

template <class T_MeshMover>
std::shared_ptr <const dolfin::Function> MeshManager<T_MeshMover>::moveMesh (const dolfin::Function& displacement, const std::string& component, const double dt)
{
    if (!displacementIsComputed_)
        computeDisplacement (displacement, component, dt);

 		// move the mesh and update mesh-dependent quantities
    oldMesh_ = * mesh_;
 		if (component != "hard")
      mesh_->move (w_);
    else
    {
/*      dolfin::Mesh tmpMesh (* mesh_);
      tmpMesh.move (w_);
      dolfin::BoundaryMesh tmpBdMesh (tmpMesh, "exterior");
      w_ = * mesh_->move (tmpBdMesh);*/
      std::vector<double> & coords (mesh_->coordinates());
      dolfin::Array<double> cs (2), xs (2);
      for (auto it (coords.begin()); it!=coords.end(); ++(++it))
      {
        xs[0] = *it;
        xs[1] = *(it+1);
        displacement.eval (cs, xs);
        *it += cs[0];
        *(it+1) += cs[1];
      }
    }
    if (component == "init")
      w_.vector()->zero();
/*    else
    {
      const std::shared_ptr<const dolfin::GenericVector> vec (w_.vector()->copy());
      w_ = dolfin::adapt (w_, mesh_);
      * w_.vector() = * vec;
    }*/
        
    //imposedDisplacement_.update ();
    displacementIsComputed_ = false;

    return dolfin::reference_to_no_delete_pointer (w_);
}

template <class T_MeshMover>
std::shared_ptr <const dolfin::Function> MeshManager<T_MeshMover>::computeDisplacement (const dolfin::Function& displacement, const std::string& component, const double dt)
{
std::cerr << "     computing displacement" << std::endl;
  if (displacementIsImposed_)
  {
std::cerr << "     imposed displacement" << std::endl;
    dolfin::Function displ (* imposedDisplacementSpace_);
    displ = imposedDisplacement_;

    displacementIsImposed_ = false;
    this->computeDisplacement (displ[0], component, dt);
    //?? togliere il [0] in qualche modo
    displacementIsImposed_ = true;

    // avoid re-computation
    displacementIsComputed_ = true;

    return dolfin::reference_to_no_delete_pointer (w_);
  }
  
  /*laplaceALE::CoefficientSpace_ux uxSp (* mesh_);
  laplaceALE::CoefficientSpace_uy uySp (* mesh_);
  dolfin::Function uxLapl (uxSp), uyLapl (uySp);
  dolfin::assign (dolfin::reference_to_no_delete_pointer (uxLapl),
                  dolfin::reference_to_no_delete_pointer (displacement[0]));
  dolfin::assign (dolfin::reference_to_no_delete_pointer (uyLapl),
                  dolfin::reference_to_no_delete_pointer (displacement[1]));
//std::cerr << wFunSp_.element()->value_rank() << ' ' << wFunSp_.str(true) << std::endl << uSp.element()->value_rank() << ' ' << uSp.str(true) << std::endl << displacement.function_space()->element()->value_rank() << ' ' << displacement.function_space()->str(true) << std::endl;
  laplacePb_->setCoefficient ("linear_form",dolfin::reference_to_no_delete_pointer(uxLapl),"ux");
  laplacePb_->setCoefficient ("linear_form",dolfin::reference_to_no_delete_pointer(uyLapl),"uy");
  */
/*  laplaceALE::CoefficientSpace u uSP (* mesh_);
  dolfin::Function uLapl (uSp);
  dolfin::assign (dolfin::reference_to_no_delete_pointer (uLapl), dolfin::reference_to_no_delete_pointer (displacement));
  laplacePb_->setCoefficient ("linear_form",dolfin::reference_to_no_delete_pointer(uLapl),"u");
*/
std::cerr << " OK fino a " << __LINE__ << std::endl;
  //dolfin::Constant nullDisplacement (0.0, 0.0);
  //aleProblem_->setCoefficient ("linear_form",dolfin::reference_to_no_delete_pointer(nullDisplacement),"u");
  //dolfin::Constant impDisplacement (0.0, 1);
  //aleProblem_->setCoefficient ("linear_form",dolfin::reference_to_no_delete_pointer(displacement),"u");
  dolfin::Function displCoeff (displCoeffSpace_);
  /*dolfin::assign (dolfin::reference_to_no_delete_pointer(displCoeff[1]), dolfin::reference_to_no_delete_pointer(displacement));
  (* displCoeff.vector()) *= 0.01;*/
std::cerr << displCoeff.function_space()->str(true) << " | " << displacement.function_space()->str(true) << std::endl;
  dolfin::assign (dolfin::reference_to_no_delete_pointer(displCoeff), dolfin::reference_to_no_delete_pointer(displacement));
    // TODO trovare un modo per passare anche la componente a cui assegnare, se serve
  aleProblem_->setCoefficient ("linear_form",dolfin::reference_to_no_delete_pointer(displCoeff),"u");
std::cerr << " OK fino a " << __LINE__ << std::endl;

  aleProblem_->solve();
  w_.vector()->zero();
  dolfin::assign (dolfin::reference_to_no_delete_pointer (w_[1]), dolfin::reference_to_no_delete_pointer (aleProblem_->solution()));
  //TODO generalize to non-vertical displacements
  * w_.vector() *= dt;

  // avoid re-computation
  displacementIsComputed_ = true;

  return dolfin::reference_to_no_delete_pointer (w_);
}

template <class T_MeshMover>
void MeshManager<T_MeshMover>::setImposedDisplacement (std::shared_ptr<const dolfin::FunctionSpace> imposedDisplacementSpace)
{
  displacementIsImposed_ = true;
  imposedDisplacementSpace_ = imposedDisplacementSpace;
}

template <class T_MeshMover>
void MeshManager<T_MeshMover>::initializeMesh (const dolfin::Mesh & mesh)
{
  //TODO
  //ma forse non si puo' fare, perche' per spostare i dof serve il FunctionSpace, che MeshManager non ha
  //motivo in piu' per farlo diventare FunctionSpaceManager invece che MeshManager
//  dolfin::dolfin_error ("MeshManager.h", "initializeMesh", "Not yet implemented (dcp)");

//    T_MeshMover::move (* this->mesh_, mesh);

  std::vector<double> newDofsCoords (wFunSp_.dofmap()->tabulate_all_coordinates (mesh));
  std::vector<double> oldDofsCoords (wFunSp_.dofmap()->tabulate_all_coordinates (* mesh_));
  std::vector<double> dofDisplacement (oldDofsCoords.size());
  std::transform (newDofsCoords.begin(), newDofsCoords.end(), oldDofsCoords.begin(), dofDisplacement.begin(), std::minus<double>());

  std::vector<double> processedDofDisplacement (dofDisplacement.size()/2);
    // NB assumiamo che in dofmap.dofs() i dof siano sempre messi a coppie, con la componente x prima e quella y dopo (quindi le coordinate si ripetono)
    // Esempio: u=(4,3) in (0.2,0.5)  ->  in dofmap.dofs() avremo 0.2,0.5,0.2,0.5  ->  dofDisplacement 0.2-newx,0.5-newy,0.2-newx,0.5-newy
  for (std::size_t k=0; k!=processedDofDisplacement.size()/2; ++k)
  {
    processedDofDisplacement [2*k] = dofDisplacement[4*k];
    processedDofDisplacement [2*k+1] = dofDisplacement[4*k+1];
  }
  dolfin::Vector displacementVec (MPI_COMM_WORLD, processedDofDisplacement.size());
  displacementVec.set_local (processedDofDisplacement);
  displacementVec.apply("insert");
  * w_.vector() = displacementVec;

	// move the mesh and update mesh-dependent quantities
  mesh_->move (w_);
  displacementIsComputed_ = false;
}

template <class T_MeshMover>
const dolfin::FunctionSpace & MeshManager<T_MeshMover>::displacementFunctionSpace () const
{
  return wFunSp_;
}

template <class T_MeshMover>
void MeshManager<T_MeshMover>::setDisplacement (const dolfin::Function & displacement)
{
  * w_.vector() = * displacement.vector();
}

template <class T_MeshMover>
void MeshManager<T_MeshMover>::print2csv (const dolfin::Function & fun, const std::string & prependToFileName, const std::string & appendToFileName, bool printAlsoDisplacement) const
{
  print2csv (fun, prependToFileName, appendToFileName, orderedDofs_);
  if (printAlsoDisplacement)
    print2csv (w_, prependToFileName, appendToFileName, orderedDisplDofs_);
}
template <class T_MeshMover>
void MeshManager<T_MeshMover>::print2csv (const dolfin::Function & fun, const std::string & prependToFileName, const std::string & appendToFileName, const std::map<std::string,std::vector<dolfin::la_index> > & dofs) const
{ 
  auto orderedDofsCoords (fun.function_space()->dofmap()->tabulate_all_coordinates (* this->mesh_));
  for (auto it (dofs.begin()); it!=dofs.end(); ++it)
  {
    std::string filename (prependToFileName+it->first+appendToFileName+".csv");
    dcp::print2csv (fun, filename, it->second.begin(), it->second.end(), orderedDofsCoords);
  }
}

} // end of namespace

#endif // MESHMANAGER_H_INCLUDE_GUARD

