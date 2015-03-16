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
#include "laplaceVec.h"

namespace geometry {

class FreeSurface : public dolfin::SubDomain
{
    bool inside (const dolfin::Array<double>& x, bool on_boundary) const
    {
        return on_boundary && 
//               !( dolfin::near (x[0], 0)
//               ||
               !(dolfin::near (x[1], 0));
//               ||
//               dolfin::near (x[0], 1));
    }
};

class FixedBoundary : public dolfin::SubDomain
{
    bool inside (const dolfin::Array<double>& x, bool on_boundary) const
    {
        return //(dolfin::near (x[0], 0) && on_boundary)
//               ||
               (dolfin::near (x[1], 0) && on_boundary);
//               ||
//               (dolfin::near (x[0], 1) && on_boundary);
    }
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

    const T_FunctionSpace& functionSpace();
    const dolfin::Mesh& mesh();

    void moveMesh (const dolfin::GenericFunction& displacement, const std::string& component="all", const double dt=1.0);
//move//    std::shared_ptr<dolfin::MeshDisplacement> moveMesh (const dolfin::GenericFunction& boundaryDisplacement);

    dolfin::Function& displacement ();
 
  protected:

    void setMeshAll (std::shared_ptr<dolfin::FunctionSpace> functionSpace, std::shared_ptr<dolfin::Mesh> mesh);

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
  problemALE_ (new dcp::LinearProblem<laplaceVec::BilinearForm, laplaceVec::LinearForm> (dolfin::reference_to_no_delete_pointer(*wFunSp_)))
/*  subElements_ (2,nullptr),
  subDofmaps_ (2,nullptr) */
  //dofmap_ (new dolfin::DofMap(std::shared_ptr<const ufc::dofmap>(new const poisson_dofmap_1), *mesh_))
  //dofmap_ (new dolfin::DofMap(*functionSpace->dofmap(),*mesh_))
  //dofmap_ (new const dolfin::DofMap(static_cast<const dolfin::DofMap&>(*functionSpace->dofmap())))
  {
    //dolfin::DofMapBuilder::build(*dofmap_,*mesh_,nullptr);
//    dofmap_ = dofmap_->create(*mesh_);
    this->moveMesh(dolfin::Constant(0,0),"initializing");
    dolfin::Constant zero (0,0);
    problemALE_->setCoefficient ("linear_form", dolfin::reference_to_no_delete_pointer (zero), "zero");
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
  problemALE_ (new dcp::LinearProblem<laplaceVec::BilinearForm, laplaceVec::LinearForm> (dolfin::reference_to_no_delete_pointer(*wFunSp_)))
/*  subElements_ (2,nullptr),
  subDofmaps_ (2,nullptr) */
  //dofmap_ (new const dolfin::DofMap(static_cast<const dolfin::DofMap&>(*functionSpace->dofmap())))
  {
    //dolfin::DofMapBuilder::build(*dofmap_,*mesh_,nullptr);
//    dofmap_ = dofmap_->create(*mesh_);
    this->moveMesh(dolfin::Constant(0,0),"initializing");
    dolfin::Constant zero (0,0);
    problemALE_->setCoefficient ("linear_form", dolfin::reference_to_no_delete_pointer (zero), "zero");
  }

template <class T_MeshMover, class T_FunctionSpace>
const T_FunctionSpace& MeshManager<T_MeshMover,T_FunctionSpace>::functionSpace()
  {
    return *functionSpace_;
  }
template <class T_MeshMover, class T_FunctionSpace>
const dolfin::Mesh& MeshManager<T_MeshMover,T_FunctionSpace>::mesh()
  {
    return *mesh_;
  }

template<class T_MeshMover, class T_FunctionSpace>
void MeshManager<T_MeshMover,T_FunctionSpace>::moveMesh (const dolfin::GenericFunction& displacement, const std::string& component, const double dt)
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
      geometry::FreeSurface freeSurface;
      geometry::FixedBoundary fixedBoundary;
      dolfin::Constant zero (0,0);
      dolfin::Constant uno (1,1);
//      geometry::ParabolicProfile parabolicProfile;
//      problemALE_->addDirichletBC (dolfin::DirichletBC (* problemALE_->functionSpace(), parabolicProfile, freeSurface));
/*      dolfin::Function displ (functionSpace_);
      displ = static_cast<const dolfin::Function&>(displacement);
      (* displ.vector()) *= dt;*/
      problemALE_->removeDirichletBC ("freeSurfaceBC");
      problemALE_->addDirichletBC (dolfin::DirichletBC (* problemALE_->functionSpace(), displacement, freeSurface, "geometric"),"freeSurfaceBC");
      problemALE_->removeDirichletBC ("fixedBoundaryBC");
      problemALE_->addDirichletBC (dolfin::DirichletBC (* problemALE_->functionSpace(), zero, fixedBoundary, "geometric"),"fixedBoundaryBC");
      problemALE_->setCoefficient ("linear_form", dolfin::reference_to_no_delete_pointer (zero), "zero");
      problemALE_->solve();
      w_ = problemALE_->solution();
      (* w_.vector()) *= dt;
dolfin::plot(w_,"la vu doppia");//dolfin::interactive();
      meshMover_.move (* mesh_, w_);
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
dolfin::interactive();
dolfin::plot(* (*(*functionSpace_)[0])[0]->mesh(), "asdjklbw [0][0] mesh dopo setMeshAll");
dolfin::interactive();
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

} // end of namespace

#endif
