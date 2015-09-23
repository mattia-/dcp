/* 
 *  Copyright (C) 2015, Ivan Fumagalli, ivanfumagalli.if@gmail.com
 * 
 *  This file is part of the DCP library
 *   
 *   The DCP library is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   The DCP library is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with the DCP library.  If not, see <http://www.gnu.org/licenses/>. 
 */ 

// HEADER file for moving-boundary geometry managing

#ifndef __GEOMETRY_H
#define __GEOMETRY_H

#include <dolfin.h>
#include <dcp/differential_problems/AbstractProblem.h>
#include <dcp/differential_problems/LinearProblem.h>
#include "laplaceVec.h"
//#include "laplaceVecSecond.h"
#include "elasticity.h"
#include "utilities.h"
#include "laplace1D.h"

//#define elasticity

namespace geometry {

class FixedBoundary : public dolfin::SubDomain
{
  public:
    bool inside (const dolfin::Array<double>& x, bool on_boundary) const
    {
        return /*(dolfin::near (x[0], 0) && on_boundary)
               ||*/
               (dolfin::near (x[1], 0) && on_boundary)
               /*||
               (dolfin::near (x[0], 1) && on_boundary)*/;
    }
};

class SlipBoundary : public dolfin::SubDomain
{
  public:
    bool inside (const dolfin::Array<double>& x, bool on_boundary) const
    {
        return ((dolfin::near (x[0], 0) && on_boundary)
                ||
                (dolfin::near (x[0], lx) && on_boundary))
//               &&
//               !dolfin::near (x[1], ly)
               &&
               !fixedBoundary_.inside(x,on_boundary);
    }
  private: 
    FixedBoundary fixedBoundary_;
};

class MovingBoundary : public dolfin::SubDomain
{
    bool inside (const dolfin::Array<double>& x, bool on_boundary) const
    {
        return on_boundary && !fixedBoundary_.inside(x,on_boundary) && (dolfin::near(x[1],ly) || !slipBoundary_.inside(x,on_boundary));
/*//               !( dolfin::near (x[0], 0)
//               ||
               !(dolfin::near (x[1], 0));
//               ||
//               dolfin::near (x[0], 1));*/
    }
    
    FixedBoundary fixedBoundary_;
    SlipBoundary slipBoundary_;
};

class ParabolicProfile : public dolfin::Expression
{
    void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x) const
    {
        values[0] = 0.5*x[0]*(1-x[0]);
        values[1] = 0.5*x[0]*(1-x[0]);
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

extern bool timeNOTcount;

class MapTgamma : public dolfin::Expression
{
// functor representing the function whose graph defines an edge of a domain

public:

  // !!! the constructor should be overridden, so that an initial value for gamma might be passed
  //void MapTgamma()

  // !!! inherited methods to be overridden

/*  void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x) const
  {
    values[0] = x[0];
    values[1] = x[1]*(1+gamma_(x[0]));
  }
*/
  void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x) const
  {
    values[0] = 0;
if (timeNOTcount)
    values[1] = x[1]*gammaDisplacement_(x[0]);// + (1.0-x[1])*(1.0-gammaDisplacement_(x[0]));
else
    values[1] = x[1]*gamma_(x[0],count);
  }

  std::size_t value_rank() const
  {
    return 1;
  }

  std::size_t value_dimension(std::size_t i) const
  {
    return 2;
  }

  // !!! new methods

//c {
  void update_count()
  {
    count++;
  }

  void init_count()
  {
    MapTgamma::count = 0;
  }
//c }

//t {
  void setTime(double t)
  {
    tOld_ = tNew_;
    tNew_ = t;
  }

  void initTime(double t)
  {
    tOld_ = t;
    tNew_ = t;
  }
//t }

/*  void set_gamma(double (*gamma)(double,int))
  {
    gamma_ = gamma;
  }
*/
  void printTimes()
  {
    std::cerr << "tOld = " << tOld_ << "  tNew = " << tNew_ << std::endl;
  }

private:

//c
  int count;
//t {
  double tNew_;
  double tOld_;
//t }
  double const PI = 3.14159265359;

/*  double gamma_(double x) const
  {
    return 0.2*sin(PI*x);
  }
*/

//c {
/*  double gamma_(double x, int param) const
  {
    return 0.2*(sin(param*PI*x)-sin((param-1)*PI*x))/(1+0.2*sin((param-1)*PI*x));
  }*/
//c {

//t {
  double gamma_ (double x, double t) const
  {
    return 1.5*t/4.0*0.4*sin(3*PI*x)+t*x;
//    return t/4.0*0.2*sin(3*PI*x);
//    return 50*x*x*(1-x)*(1-x)*(x-0.3)*(x-0.7)*(t-1)*(t-3)*(1+t/2);
    //return 75*x*x*(1-x)*(1-x)*(x-0.3)*(x-0.7)*t*(2-t)*(4-t)/(1+t);
    if (t==0)
      return 0;
    return 200*(-(x-t-0.5)*(x-t-0.5)*(x-t-1.0)*(x-t-1.0))*(x>=t+0.5 ? 1 : 0)*(x<=t+1.0 ? 1 : 0);
  }

  double gammaDisplacement_ (double x) const
  {
    return (gamma_(x,tNew_)-gamma_(x,tOld_))/(1+gamma_(x,tOld_));
  }
//t }

}; // end of class gamma

/*template < T_FunSp >
class FunctionSpace : public dolfin::FunctionSpace {

public:

  FunctionSpace(T_FunSp & funSp) : funSp_(funSp) {};

  /// Return mesh
  ///
  /// *Returns*
  ///     _Mesh_
  ///         The mesh.
  std::shared_ptr<const Mesh> mesh() const
  {
    return mesh_;
  }

  /// Return finite element
  ///
  /// *Returns*
  ///     _FiniteElement_
  ///         The finite element.
  std::shared_ptr<const FiniteElement> element() const
  {
    return element_;
  }

  /// Return dofmap
  ///
  /// *Returns*
  ///     _GenericDofMap_
  ///         The dofmap.
  std::shared_ptr<const GenericDofMap> dofmap() const
  {
    return dofmap_;
  }

protected:

  poisson::FunctionSpace funSp_;

private:
  // All the following members are private in the base class
  // dolfin::FunctionSpace, then they have to be defined from scratch
  
  // The mesh
  std::shared_ptr<const Mesh> mesh_;

  // The finite element
  std::shared_ptr<const FiniteElement> element_;

  // The dofmap
  std::shared_ptr<const GenericDofMap> dofmap_;

  // The component (for subspaces)
  std::vector<std::size_t> component_;

  // Cache of subspaces
  mutable std::map<std::vector<std::size_t>,
    std::shared_ptr<FunctionSpace> > subspaces_;

}; // end of class FunctionSpace
*/

template <class T_MeshMover, class T_FunctionSpace>
class MeshManager
{

  public:

    typedef T_MeshMover MeshMover;
    typedef T_FunctionSpace FunctionSpace;

    MeshManager () = delete;
  
    MeshManager (const std::shared_ptr<T_MeshMover> meshMover, std::shared_ptr<T_FunctionSpace> functionSpace);
    MeshManager (const T_MeshMover& meshMover, T_FunctionSpace& functionSpace);
    // !!! servono altri costruttori?

    const T_FunctionSpace& functionSpace() const;
    std::shared_ptr<const dolfin::Mesh> mesh() const;

    void moveMesh (const dolfin::Function& displacement, const std::string& component="all", const double dt=1.0);
    void moveMeshOld (const dolfin::Function& displacement, const std::string& component="all", const double dt=1.0);
//move//    std::shared_ptr<dolfin::MeshDisplacement> moveMesh (const dolfin::GenericFunction& boundaryDisplacement);

    dolfin::Function& displacement ();
 
  protected:

    void setMeshAll (std::shared_ptr<dolfin::FunctionSpace> functionSpace, std::shared_ptr<dolfin::Mesh> mesh);
    void smooth (std::shared_ptr<dolfin::Mesh> mesh);
    dolfin::Function preSmooth (const dolfin::Function& displacement, const std::size_t component=1);

    T_MeshMover meshMover_;
    std::shared_ptr<T_FunctionSpace> functionSpace_;
    std::shared_ptr<dolfin::Mesh> mesh_;
    std::shared_ptr<const dolfin::FiniteElement> element_;
    std::shared_ptr<const dolfin::GenericDofMap> dofmap_;
    std::shared_ptr<dolfin::FunctionSpace> wFunSp_;
    dolfin::Function w_;
    std::shared_ptr<dcp::AbstractProblem> problemALE_;
/*    std::vector<std::shared_ptr<dolfin::FiniteElement> > subElements_;
    std::vector<std::shared_ptr<dolfin::GenericDofMap> > subDofmaps_;
*/
  dolfin::FacetFunction<std::size_t> meshFacets_;
    dolfin::MeshFunction<bool> on_boundary_;
    dolfin::MeshFunction<bool> on_to_be_neglected_boundary_;
    ImposedDisplacement imposedDisplacement_;
    //ImposedFromFile imposedDisplacement_;
    dolfin::Mesh intervalMesh_;
    laplace1D::FunctionSpace intervalFunSp_;
    dolfin::Function intervalSol_;
    FunctionToExpression preSmoothed_;

};

// IMPLEMENTATION

template <class T_MeshMover, class T_FunctionSpace>
MeshManager<T_MeshMover,T_FunctionSpace>::MeshManager(const std::shared_ptr<T_MeshMover> meshMover, std::shared_ptr<T_FunctionSpace> functionSpace) :
  meshMover_ (*meshMover),
  functionSpace_ (functionSpace),
  mesh_ (new dolfin::Mesh(*functionSpace->mesh())),
  element_ (new dolfin::FiniteElement(*functionSpace->element())),
  dofmap_ (functionSpace->dofmap()->copy()),
//  problemALE_ (dolfin::reference_to_no_delete_pointer(*mesh_), dolfin::reference_to_no_delete_pointer(laplaceVec::FunctionSpace(*mesh_))),
//  w_ (* problemALE_->functionSpace())
  wFunSp_ (new laplaceVec::FunctionSpace(*mesh_)),
  w_ (*wFunSp_),
 // problemALE_ (new dcp::LinearProblem<laplaceVec::BilinearForm, laplaceVec::LinearForm> (dolfin::reference_to_no_delete_pointer(*mesh_), dolfin::reference_to_no_delete_pointer(*wFunSp_)))
  problemALE_ (new dcp::LinearProblem<laplaceVec::BilinearForm, laplaceVec::LinearForm> (dolfin::reference_to_no_delete_pointer(*wFunSp_))),
/*  subElements_ (2,nullptr),
  subDofmaps_ (2,nullptr) */
  //dofmap_ (new dolfin::DofMap(std::shared_ptr<const ufc::dofmap>(new const poisson_dofmap_1), *mesh_))
  //dofmap_ (new dolfin::DofMap(*functionSpace->dofmap(),*mesh_))
  //dofmap_ (new const dolfin::DofMap(static_cast<const dolfin::DofMap&>(*functionSpace->dofmap())))
  meshFacets_ (* mesh_),
  on_boundary_ (* mesh_,0),
  on_to_be_neglected_boundary_ (*mesh_, 0),
  intervalMesh_ (dolfin::UnitIntervalMesh(ny)),
  intervalFunSp_ (intervalMesh_),
  intervalSol_ (intervalFunSp_),//laplace1D::FunctionSpace(dolfin::UnitIntervalMesh(ny))),
  preSmoothed_ (dolfin::reference_to_no_delete_pointer(intervalSol_),2,1)
  {
    //dolfin::DofMapBuilder::build(*dofmap_,*mesh_,nullptr);
//    dofmap_ = dofmap_->create(*mesh_);
    dolfin::Function initFun (* wFunSp_);
    initFun = dolfin::Constant(0,0);
//    this->moveMesh(initFun,"initializing");
    dolfin::Constant zero (0,0);
    problemALE_->setCoefficient ("linear_form", dolfin::reference_to_no_delete_pointer (zero), "zero");
  meshFacets_.set_all (0);
  TopBoundary freeSurface;
  LateralWall lateralWall;
  BottomBoundary bottomBoundary;
  freeSurface.mark (meshFacets_, 1);
  lateralWall.mark (meshFacets_, 2);
  bottomBoundary.mark (meshFacets_, 3);

dolfin::plot(*mesh_, "pre smoothing");//dolfin::interactive();
    // Make sure we have cell-facet connectivity
    mesh_->init(mesh_->topology().dim(), mesh_->topology().dim() - 1);
    // Make sure we have vertex-edge connectivity
    mesh_->init(0, 1);
    // Make sure the mesh is ordered
    mesh_->order();

    dolfin::BoundaryMesh boundary (*mesh_, "exterior");
/*    std::vector<double>& coords (boundary.coordinates());
    const std::vector<unsigned int>& cells (boundary.cells());*/

    // Mark vertices on the boundary so we may skip the internal ones
    const dolfin::MeshFunction<std::size_t> vertex_map = boundary.entity_map(0);
    on_boundary_  = false;
    on_to_be_neglected_boundary_  = false;
    if (boundary.num_vertices() > 0)
    {
      for (dolfin::VertexIterator v(boundary); !v.end(); ++v)
      {
          on_boundary_ [vertex_map[*v]] = lateralWall.inside(dolfin::Array<double>(2,const_cast<double*>(v->x())),true);
//          on_to_be_neglected_boundary_ [vertex_map[*v]] = ! lateralWall.inside(dolfin::Array<double>(2,const_cast<double*>(v->x())),on_boundary_ [vertex_map[*v]]);
          on_to_be_neglected_boundary_ [vertex_map[*v]] = lateralWall.inside(dolfin::Array<double>(2,const_cast<double*>(v->x())),on_boundary_ [vertex_map[*v]])
                                                        &&
                                                        ( freeSurface.inside(dolfin::Array<double>(2,const_cast<double*>(v->x())),on_boundary_ [vertex_map[*v]])
                                                          || 
                                                          bottomBoundary.inside(dolfin::Array<double>(2,const_cast<double*>(v->x())),on_boundary_ [vertex_map[*v]])
                                                        );
      }
    }
/*freeSurface.mark (on_boundary_ ,false);
bottomBoundary.mark (on_boundary_ ,false);*/

  }

template <class T_MeshMover, class T_FunctionSpace>
MeshManager<T_MeshMover,T_FunctionSpace>::MeshManager (const T_MeshMover& meshMover, T_FunctionSpace& functionSpace) :
  meshMover_ (meshMover),
  functionSpace_ (&functionSpace),
  mesh_ (new dolfin::Mesh(*functionSpace.mesh())),
  element_ (new const dolfin::FiniteElement(*functionSpace->element())),
  dofmap_ (functionSpace->dofmap()->copy()),
//  problemALE_ (dolfin::reference_to_no_delete_pointer(*mesh_), dolfin::reference_to_no_delete_pointer(laplaceVec::FunctionSpace(*mesh_))),
//  w_ (* problemALE_->functionSpace())
  wFunSp_ (new laplaceVec::FunctionSpace(*mesh_)),
  w_ (*wFunSp_),
 // problemALE_ (new dcp::LinearProblem<laplaceVec::BilinearForm, laplaceVec::LinearForm> (dolfin::reference_to_no_delete_pointer(*mesh_), dolfin::reference_to_no_delete_pointer(*wFunSp_)))
  problemALE_ (new dcp::LinearProblem<laplaceVec::BilinearForm, laplaceVec::LinearForm> (dolfin::reference_to_no_delete_pointer(*wFunSp_))),
/*  subElements_ (2,nullptr),
  subDofmaps_ (2,nullptr) */
  //dofmap_ (new const dolfin::DofMap(static_cast<const dolfin::DofMap&>(*functionSpace->dofmap())))
  meshFacets_ (* mesh_)
  {
//TODO: sistemare: per ora uso l'altro e bon
std::cerr << "!!! USA L'ALTRO COSTRUTTORE DI MeshManager" << std::endl;
exit(1);
    //dolfin::DofMapBuilder::build(*dofmap_,*mesh_,nullptr);
//    dofmap_ = dofmap_->create(*mesh_);
    dolfin::Function initFun (* wFunSp_);
    initFun = dolfin::Constant(0,0);
//    this->moveMesh(initFun,"initializing");
    dolfin::Constant zero (0,0);
    problemALE_->setCoefficient ("linear_form", dolfin::reference_to_no_delete_pointer (zero), "zero");
  meshFacets_.set_all (0);
  TopBoundary freeSurface;
  LateralWall lateralWall;
  BottomBoundary bottomBoundary;
  freeSurface.mark (meshFacets_, 1);
  lateralWall.mark (meshFacets_, 2);
  bottomBoundary.mark (meshFacets_, 3);
  }

template <class T_MeshMover, class T_FunctionSpace>
const T_FunctionSpace& MeshManager<T_MeshMover,T_FunctionSpace>::functionSpace() const
  {
    return *functionSpace_;
  }
template <class T_MeshMover, class T_FunctionSpace>
std::shared_ptr<const dolfin::Mesh> MeshManager<T_MeshMover,T_FunctionSpace>::mesh() const
  {
    return mesh_;
  }

template<class T_MeshMover, class T_FunctionSpace>
void MeshManager<T_MeshMover,T_FunctionSpace>::moveMeshOld (const dolfin::Function& displacement, const std::string& component, const double dt)
{
#ifdef elasticity
    wFunSp_.reset (new elasticity::FunctionSpace(*mesh_));
#else
    wFunSp_.reset (new laplaceVec::FunctionSpace(*mesh_));
#endif
    dolfin::Function w (displacement);//*wFunSp_);
//    (* w.vector()) = (* displacement.vector());
    (* w.vector()) *= dt;
    
/*    mesh_->move (w);
    mesh_->smooth_boundary (5);
//    mesh_->smooth (5);*/
    this->smooth (mesh_);
dolfin::plot(*mesh_, "moved mesh"); //dolfin::interactive();

    setMeshAll(functionSpace_, mesh_);
    w_ = w;
}

template<class T_MeshMover, class T_FunctionSpace>
void MeshManager<T_MeshMover,T_FunctionSpace>::moveMesh (const dolfin::Function& displacement, const std::string& component, const double dt)
{
  //  std::cout << "muovo la mesh di " << functionSpace_->str(true) << std::endl;
  //  std::cout << element_->value_rank() << " " << element_->value_dimension(1) << " " << element_->space_dimension() << " " << element_->num_sub_elements() << " " << element_->signature() << std::endl;
  if (component == "all" || component == "initializing")
  {
//      meshMover_.move (const_cast<dolfin::Mesh&>(* functionSpace_->mesh()), displacement);
      meshMover_.move (* mesh_, displacement);
  }
  else if (component == "normal")
  {
      geometry::MovingBoundary movingBoundary;
      geometry::FixedBoundary fixedBoundary;
      geometry::SlipBoundary slipBoundary;
      dolfin::Constant zero (0,0);
      dolfin::Constant uno (1,1);
      dolfin::Constant zeroComp (0);
      dolfin::Constant penalty (1.6e3);
      dolfin::Constant volMin (3.14159265 * mesh_->rmin()*mesh_->rmin());
      dolfin::Constant volMax (0.5 * 3.14159265 * mesh_->hmax()*mesh_->hmax()*0.25);
//      geometry::ParabolicProfile parabolicProfile;
//      problemALE_->addDirichletBC (dolfin::DirichletBC (* problemALE_->functionSpace(), parabolicProfile, freeSurface));
/*      dolfin::Function displ (functionSpace_);
      displ = static_cast<const dolfin::Function&>(displacement);
      (* displ.vector()) *= dt;*/
TopBoundary freeSurface;
/*  dolfin::FacetFunction<std::size_t> meshFacets (*mesh_);
  meshFacets.set_values (std::vector<std::size_t>(meshFacets_.values(),meshFacets_.values()+meshFacets_.size()));
//for (std::size_t i=0; i!= meshFacets.size(); ++i) std::cerr << meshFacets.values()[i] << ", "; std::cerr << std::endl;
  wFunSp_.reset (new laplaceVec::FunctionSpace(*mesh_));
  problemALE_.reset (new dcp::LinearProblem<laplaceVec::BilinearForm, laplaceVec::LinearForm> (dolfin::reference_to_no_delete_pointer(*wFunSp_)));
      problemALE_->removeDirichletBC ("movingBoundaryBC");
//      problemALE_->addDirichletBC (dolfin::DirichletBC (* problemALE_->functionSpace(), uno, meshFacets_, 1),"movingBoundaryBC");
      problemALE_->addDirichletBC (dolfin::DirichletBC (* problemALE_->functionSpace(), uno, freeSurface),"movingBoundaryBC");
//      problemALE_->addDirichletBC (dolfin::DirichletBC (* problemALE_->functionSpace(), uno, meshFacets_, 1),"movingBoundaryBC");
dolfin::plot(meshFacets,"meshFacets per BC");
dolfin::plot(meshFacets_,"meshFacets_ per BC");
//      problemALE_->addDirichletBC (dolfin::DirichletBC (* (* problemALE_->functionSpace())[1], displacement[1], movingBoundary, "geometric"),"movingBoundaryBC");
      problemALE_->removeDirichletBC ("fixedBoundaryBC");
      problemALE_->addDirichletBC (dolfin::DirichletBC (* problemALE_->functionSpace(), zero, fixedBoundary, "geometric"),"fixedBoundaryBC");
//      problemALE_->addDirichletBC (dolfin::DirichletBC (* problemALE_->functionSpace(), displacement, fixedBoundary, "geometric"),"fixedBoundaryBC");
      problemALE_->removeDirichletBC ("slipBoundaryBC");
      problemALE_->addDirichletBC (dolfin::DirichletBC (* (* problemALE_->functionSpace())[0], zeroComp, slipBoundary, "geometric"),"slipBoundaryBC");
//      problemALE_->addDirichletBC (dolfin::DirichletBC (* problemALE_->functionSpace(), displacement, slipBoundary, "geometric"),"slipBoundaryBC");
      problemALE_->setCoefficient ("linear_form", dolfin::reference_to_no_delete_pointer (zero), "zero");
      problemALE_->solve();
      w_ = problemALE_->solution();
      (* w_.vector()) *= dt;
      meshMover_.move (* mesh_, w_);
dolfin::plot(w_,"la vu doppia");//dolfin::interactive();
dolfin::plot(* wFunSp_->mesh(), "wFunSp_ mesh");dolfin::interactive();*/
std::cerr << "OK fino a " << __LINE__ << std::endl;
#ifdef elasticity
  wFunSp_.reset (new elasticity::FunctionSpace(*mesh_));
  problemALE_.reset (new dcp::LinearProblem<elasticity::BilinearForm, elasticity::LinearForm> (dolfin::reference_to_no_delete_pointer(*wFunSp_)));
#else
  wFunSp_.reset (new laplaceVec::FunctionSpace(*mesh_));
  problemALE_.reset (new dcp::LinearProblem<laplaceVec::BilinearForm, laplaceVec::LinearForm> (dolfin::reference_to_no_delete_pointer(*wFunSp_)));
#endif
  problemALE_->setIntegrationSubdomain ("bilinear_form",
        dolfin::reference_to_no_delete_pointer (meshFacets_), dcp::SubdomainType::BOUNDARY_FACETS);
  problemALE_->setIntegrationSubdomain ("linear_form",
        dolfin::reference_to_no_delete_pointer (meshFacets_), dcp::SubdomainType::BOUNDARY_FACETS);
std::cerr << "OK fino a " << __LINE__ << std::endl;
#ifdef elasticity
  problemALE_->setCoefficient ("bilinear_form", dolfin::reference_to_no_delete_pointer (penalty), "penalty");
  problemALE_->setCoefficient ("linear_form", dolfin::reference_to_no_delete_pointer (penalty), "penalty");
#endif
/*  problemALE_->setCoefficient ("bilinear_form", dolfin::reference_to_no_delete_pointer (volMin), "volMin");
  problemALE_->setCoefficient ("bilinear_form", dolfin::reference_to_no_delete_pointer (volMax), "volMax");
  problemALE_->setCoefficient ("linear_form", dolfin::reference_to_no_delete_pointer (volMin), "volMin");
  problemALE_->setCoefficient ("linear_form", dolfin::reference_to_no_delete_pointer (volMax), "volMax");*/
std::cerr << "OK fino a " << __LINE__ << std::endl;
std::unordered_map<std::size_t, double> bdval;
      problemALE_->removeDirichletBC ("fixedBoundaryBC");
//      problemALE_->addDirichletBC (dolfin::DirichletBC (* problemALE_->functionSpace(), zero, fixedBoundary, "geometric"),"fixedBoundaryBC");
      dolfin::DirichletBC dirBCfixed (* problemALE_->functionSpace(), zero, this->meshFacets_, 3, "geometric");
dirBCfixed.get_boundary_values(bdval);
      problemALE_->addDirichletBC (dirBCfixed,"fixedBoundaryBC");
//      problemALE_->addDirichletBC (dolfin::DirichletBC (* problemALE_->functionSpace(), displacement, fixedBoundary, "geometric"),"fixedBoundaryBC");
      problemALE_->removeDirichletBC ("slipBoundaryBC");
//      problemALE_->addDirichletBC (dolfin::DirichletBC (* (* problemALE_->functionSpace())[0], zeroComp, slipBoundary, "geometric"),"slipBoundaryBC");
      dolfin::DirichletBC dirBCslip (* (* problemALE_->functionSpace())[0], zeroComp, this->meshFacets_, 2, "geometric");
dirBCfixed.get_boundary_values(bdval);
      problemALE_->addDirichletBC (dirBCslip,"slipBoundaryBC");
//      problemALE_->addDirichletBC (dolfin::DirichletBC (* problemALE_->functionSpace(), displacement, slipBoundary, "geometric"),"slipBoundaryBC");
      problemALE_->removeDirichletBC ("FSBoundaryBC");
      problemALE_->removeDirichletBC ("WallBoundaryBC");
/*dolfin::plot(this->meshFacets_);
std::vector<std::size_t> facets;
for (dolfin::FacetIterator facet(* problemALE_->functionSpace()->mesh()); !facet.end(); ++facet)
  {
    if (meshFacets_[*facet] == 0)
      facets.push_back(facet->index());
  }
std::cerr << "# facets = " << facets.size() << std::endl;
dolfin::Constant prova(1,1);*/
/*this->imposedDisplacement_.setMesh (* mesh_);
this->imposedDisplacement_.loadFun();*/
      dolfin::Function displ (* wFunSp_);
std::cerr << " ciao " << displ.value_size() << std::endl;
//      displ = this->imposedDisplacement_;
//dolfin::Function preSmoothedFun (preSmooth (displacement));
//std::cerr << " ciao " << preSmoothedFun.value_rank() << std::endl;
//FunctionToExpression preSmoothed (preSmoothedFun,2,1);
//std::cerr << (preSmooth (displacement))(0.000001,0);
//FunctionToExpression preSmoothed (preSmooth (displacement),displacement.value_size(),1);
//preSmooth (displacement);
std::cerr << " ciao " << preSmoothed_.value_size() << std::endl;
displ = preSmoothed_;
std::cerr << " ciao " << displ.value_size() << std::endl;
dolfin::plot (preSmoothed_, *mesh_, "preSmoothed"); //dolfin::interactive();
dolfin::plot (displ, "preSmoothed displ"); //dolfin::interactive();

//      dolfin::DirichletBC dirBC (* problemALE_->functionSpace(), this->imposedDisplacement_, this->meshFacets_,0, "geometric");
//      dolfin::DirichletBC dirBCfs (* (* problemALE_->functionSpace())[1], this->imposedDisplacement_, this->meshFacets_,1, "topological");
      dolfin::DirichletBC dirBCfs (* (* problemALE_->functionSpace())[1], displacement[1], this->meshFacets_,1, "topological");
   //   dolfin::DirichletBC dirBCwall (* problemALE_->functionSpace(), displ, this->meshFacets_,2, "geometric");
      dolfin::DirichletBC dirBCwall (* (* problemALE_->functionSpace())[0], zeroComp, this->meshFacets_,2, "geometric");
  //    dolfin::DirichletBC dirBCfs (* problemALE_->functionSpace(), displacement, this->meshFacets_,1, "topological");
dirBCfs.get_boundary_values(bdval);
dirBCwall.get_boundary_values(bdval);
  //???!!! se tolgo questa riga sopra, pare che la DirichletBC non inizializzi correttamente il suo membro _facets...  ,:-[
/*std::cerr << "# facets dentro la BC = " << bdval.size();
std::cerr << "# facets dentro la BC = " << dirBC.markers().size() << std::endl;*/
      problemALE_->addDirichletBC (dirBCfs,"FSBoundaryBC");
      problemALE_->addDirichletBC (dirBCwall,"WallBoundaryBC");
/*problemALE_->dirichletBC("FSBoundaryBC").get_boundary_values(bdval);
std::cerr << "# facets dentro la BC dentro il problem = " << bdval.size();
for(auto it=bdval.begin(); it!=bdval.end(); it++) std::cerr << it->first << " " << it->second << std::endl;
std::cerr << "# facets dentro la BC dentro il problem = " << problemALE_->dirichletBC("FSBoundaryBC").markers().size() << std::endl;*/
std::cerr << "OK fino a " << __LINE__ << std::endl;
      problemALE_->setCoefficient ("linear_form", dolfin::reference_to_no_delete_pointer (zero), "zero");
std::cerr << "OK fino a " << __LINE__ << std::endl;
#ifdef elasticity
      problemALE_->setCoefficient ("linear_form", dolfin::reference_to_no_delete_pointer (displacement), "datum");
#endif
//      problemALE_->setCoefficient ("bilinear_form", dolfin::reference_to_no_delete_pointer (*(new ALEStiffnessVec)),"stiffnessVec");
//      problemALE_->setCoefficient ("linear_form", dolfin::reference_to_no_delete_pointer (*(new ALEStiffnessVec)),"stiffnessVec");
std::cerr << "OK fino a " << __LINE__ << std::endl;
      problemALE_->solve();
      w_ = problemALE_->solution();
      (* w_.vector()) *= dt;
      meshMover_.move (* mesh_, w_);
//this->imposedDisplacement_.update();
this->preSmoothed_.updateRef(dt);
//??? PERCHE' CAVOLO NON VA CON 'STA RIGA??

//mesh_->smooth_boundary(10,false);
//mesh_->smooth(5);

dolfin::plot(w_,"la vu doppia second");//dolfin::interactive();
dolfin::plot(* wFunSp_->mesh(), "wFunSp_ mesh second");//dolfin::interactive();
  }
//  else if (component == "tangential")
  else
  {
      std::cerr << "Unadmissible component choice (options are: all, normal, tangential)"
                << "\t in void MeshManager<T_MeshMover,T_FunctionSpace>::moveMesh (const dolfin::GenericFunction& displacement, const std::string& component)"
                << std::endl;
      exit(10);
  }
  //  *functionSpace_ = T_FunctionSpace(mesh_,element_,dofmap_);
  /*  for (int i=0; i!=2; i++)
  {
  //    std::cerr << "componente " << *iter << "  " << comp->str(false) << std::endl;
  std::cerr << "componente " << i << "  " << (*functionSpace_)[i]->str(true) << std::endl;
  std::string title("mesh della componente "); title += std::to_string(i);
  dolfin::plot(*(*functionSpace_)[i]->mesh(),title.c_str());
  dolfin::interactive();
  // sistema le componenti
  subElements_[i] = std::shared_ptr<dolfin::FiniteElement> (new dolfin::FiniteElement(*(*functionSpace_)[i]->element()));
  subDofmaps_[i]  = std::shared_ptr<dolfin::GenericDofMap>((*functionSpace_)[i]->dofmap()->copy());
  std::cout << subElements_[i]->value_rank() << " " << subElements_[i]->value_dimension(1) << " " << subElements_[i]->space_dimension() << " " << subElements_[i]->num_sub_elements() << " " << subElements_[i]->signature() << std::endl;
  (*(*functionSpace_)[i]) = T_FunctionSpace(mesh_,subElements_[i],subDofmaps_[i]);
  title += " dopo";
  dolfin::plot(*(*functionSpace_)[i]->mesh(),title.c_str());
  dolfin::interactive();
  }
  */
//dolfin::plot(*functionSpace_->mesh(), "asdjklbw mesh");
//dolfin::interactive();
dolfin::plot(* (*(*functionSpace_)[0])[0]->mesh(), "asdjklbw [0][0] mesh");
//dolfin::interactive();
  setMeshAll(functionSpace_, mesh_);
  setMeshAll(problemALE_->functionSpace(), mesh_);
if (component == "initializing") {
dolfin::plot(*functionSpace_->mesh(), "asdjklbw mesh dopo setMeshAll");
//dolfin::interactive();
dolfin::plot(* (*(*functionSpace_)[0])[0]->mesh(), "asdjklbw [0][0] mesh dopo setMeshAll");
//dolfin::interactive();
/*  * functionSpace_ = dolfin::FunctionSpace (mesh_,functionSpace_->element(),functionSpace_->dofmap()->create(*mesh_));
dolfin::plot(*functionSpace_->mesh(), "asdjklbw mesh dopo ricostruzione");
dolfin::interactive();
dolfin::plot(* (*(*functionSpace_)[0])[0]->mesh(), "asdjklbw [0][0] mesh dopo ricostruzione");
dolfin::interactive();*/
}
  // TODO forse questo casino dei subspace e' eliminabile, se T_FunctionSpace potesse essere diverso da dolfin::FunctionSpace,
  // pero' servirebbe che AbstractProblem fosse templatizzato sul tipo di functionSpace_, cosi' da poter creare un getter che
  // restituisca un puntatore / una referenza (non const) al tipo "vero" del function space
  // (non ne sono certo, ma  *functionSpace_=T_FunctionSpace(mesh_,element_,dofmap_);  dovrebbe funzionare)
}

//move//  template<class T_MeshMover, class T_FunctionSpace>
//move//  std::shared_ptr<dolfin::MeshDisplacement> MeshManager<T_MeshMover,T_FunctionSpace>::moveMesh (const dolfin::GenericFunction& boundaryDisplacement)
//move//  {
//move//    dolfin::Mesh tmpMesh(*mesh_);
//move//    meshMover_.move(tmpMesh,boundaryDisplacement);
//move//  //  tmpMesh.move(boundaryDisplacement);
//move//  dolfin::plot(tmpMesh,"tmpMesh");
//move//  dolfin::interactive();
//move//  
//move//    std::shared_ptr<dolfin::MeshDisplacement> displacement (new dolfin::MeshDisplacement(*(meshMover_.move(*mesh_,tmpMesh))));
//move//  dolfin::plot(*mesh_,"mesh_");
//move//  dolfin::interactive();
//move//    setMeshAll(functionSpace_,mesh_);
//move//  
//move//    return displacement;
//move//  }
 
template<class T_MeshMover, class T_FunctionSpace>
dolfin::Function& MeshManager<T_MeshMover,T_FunctionSpace>::displacement ()
{
  return w_;
}

template<class T_MeshMover, class T_FunctionSpace>
void MeshManager<T_MeshMover,T_FunctionSpace>::setMeshAll (std::shared_ptr<dolfin::FunctionSpace> functionSpace, std::shared_ptr<dolfin::Mesh> mesh)
{

  for (std::size_t i=0; i!=functionSpace->element()->num_sub_elements(); i++)
  {
    setMeshAll ((*functionSpace)[i], mesh);
  }
  *functionSpace = dolfin::FunctionSpace(mesh,functionSpace->element(),functionSpace->dofmap());
//  std::shared_ptr<dolfin::Mesh> pMesh (std::const_pointer_cast<dolfin::Mesh> (functionSpace->mesh()));
//  *pMesh = *mesh;

//  functionSpace.reset (new T_FunctionSpace(*mesh));

}

template<class T_MeshMover, class T_FunctionSpace>
void MeshManager<T_MeshMover,T_FunctionSpace>::smooth (std::shared_ptr<dolfin::Mesh> mesh)
{
    const std::size_t d = mesh->geometry().dim();
    std::vector<double> xx(d);
std::size_t num_iterations (3);
    for (std::size_t iteration = 0; iteration < num_iterations; iteration++)
    {
      for (dolfin::VertexIterator v(*mesh); !v.end(); ++v)
      {
        // Skip internal vertices
        if ( (! on_boundary_ [*v]) || (on_to_be_neglected_boundary_ [*v]) )
          continue;

        // Get coordinates of vertex
        double* x = mesh->geometry().x(v->index());
        const dolfin::Point p = v->point();
std::cerr << "Vertice corrente: " << x[0] << " , " << x[1] << std::endl;

        // Compute center of mass of neighboring vertices
        for (std::size_t i = 0; i < d; i++) xx[i] = 0.0;
          std::size_t num_neighbors = 0;
        for (dolfin::EdgeIterator e(*v); !e.end(); ++e)
        {
          // Get the other vertex
          dolfin_assert(e->num_entities(0) == 2);
          std::size_t other_index = e->entities(0)[0];
          if (other_index == v->index())
            other_index = e->entities(0)[1];
          // Create the vertex
          dolfin::Vertex vn(*mesh, other_index);
          // Skip the vertex itself
          if (v->index() == vn.index())
            continue;
          // Skip internal vertices
          if ( ! on_boundary_ [vn] )
            continue;
          num_neighbors += 1;
          // Compute center of mass
          const double* xn = vn.x();
          for (std::size_t i = 0; i < d; i++)
            xx[i] += xn[i];
std::cerr << "    altro vertice: " << xn[0] << " , " << xn[1] << std::endl;
        }
        for (std::size_t i = 0; i < d; i++)
          xx[i] /= static_cast<double>(num_neighbors);
          
        // Move vertex
        for (std::size_t i = 0; i < d; i++)
          x[i] = xx[i];
      }
dolfin::plot(*mesh, "post smoothing");d//olfin::interactive();
    }
    mesh->smooth(5);
dolfin::plot(*mesh, "post dolfin smoothing");//dolfin::interactive();
}

template<class T_MeshMover, class T_FunctionSpace>
dolfin::Function MeshManager<T_MeshMover,T_FunctionSpace>::preSmooth (const dolfin::Function& displacement, const std::size_t component)
{
/*    dolfin::IntervalMesh intervalMesh (ny,0,1);
    laplace1D::FunctionSpace intervalFunSpace (intervalMesh);*/
    laplace1D::BilinearForm a (* intervalSol_.function_space(),* intervalSol_.function_space());
    laplace1D::LinearForm F (* intervalSol_.function_space());
    dolfin::Constant zero (0);
    F.set_coefficient ("zero",dolfin::reference_to_no_delete_pointer(zero));

    FirstExtr firstExtr;
    SecondExtr secondExtr;
    dolfin::DirichletBC bc1 (* intervalSol_.function_space(),zero,firstExtr);//displacement[component],firstExtr);
std::cerr << " OK fino a " << __LINE__ << std::endl;
    FunctionToExpressionBis displ(dolfin::reference_to_no_delete_pointer(displacement),component,preSmoothed_.length());
dolfin::plot(displ,* intervalSol_.function_space()->mesh(),"displ in preSmooth");//dolfin::interactive();
    dolfin::DirichletBC bc2 (* intervalSol_.function_space(),displ,secondExtr);
    std::vector<const dolfin::DirichletBC*> bcs ({& bc2, & bc1});
        // l'ordine conta in questo vector!
        
    dolfin::solve (a==F,intervalSol_,bcs);
dolfin::plot(intervalSol_,"intervalSol");//dolfin::interactive();

std::cerr << "salut" << intervalSol_.value_size() << std::endl;
//    FunctionToExpression returnSol (intervalSol,displacement.value_size(),component);

    return intervalSol_;//returnSol;
}

} // end of namespace

#endif
