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
#include "poisson.h"

namespace geometry {

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
    values[1] = x[1]*gamma_(x[0],count);
//t    values[1] = x[1]*gammaDisplacement_(x[0]);
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

  void update_count()
  {
    count++;
  }

  void init_count()
  {
    MapTgamma::count = 0;
  }

/*//t  void setTime(double t)
  {
    tOld_ = tNew_;
    tNew_ = t;
  }

  void initTime(double t)
  {
    tOld_ = t;
    tNew_ = t;
  }
*///t

/*  void set_gamma(double (*gamma)(double,int))
  {
    gamma_ = gamma;
  }
*/

private:

  int count;
//t  double tNew_;
//t  double tOld_;
  double const PI = 3.14159265359;

  double gamma_(double x) const
  {
    return 0.2*sin(PI*x);
  }

  double gamma_(double x, int param) const
  {
    return 0.2*(sin(param*PI*x)-sin((param-1)*PI*x))/(1+0.2*sin((param-1)*PI*x));
  }

/*//t  double gamma_ (double x, double t) const
  {
    return t/4.0*0.2*sin(3*PI*x);
  }

  double gammaDisplacement_ (double x) const
  {
    return gamma_(x,tNew_)-gamma_(x,tOld_);
  }
*///t

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

    void moveMesh (const dolfin::GenericFunction& displacement);
 
  private:

    T_MeshMover meshMover_;
    std::shared_ptr<T_FunctionSpace> functionSpace_;
    std::shared_ptr<dolfin::Mesh> mesh_;
    std::shared_ptr<const dolfin::FiniteElement> element_;
    std::shared_ptr<const dolfin::GenericDofMap> dofmap_;
    std::vector<std::shared_ptr<dolfin::FiniteElement> > subElements_;
    std::vector<std::shared_ptr<dolfin::GenericDofMap> > subDofmaps_;

};

// IMPLEMENTATION

template <class T_MeshMover, class T_FunctionSpace>
MeshManager<T_MeshMover,T_FunctionSpace>::MeshManager(const std::shared_ptr<T_MeshMover> meshMover, std::shared_ptr<T_FunctionSpace> functionSpace) :
  meshMover_ (*meshMover),
  functionSpace_ (functionSpace),
  mesh_ (new dolfin::Mesh(*functionSpace->mesh())),
  element_ (new dolfin::FiniteElement(*functionSpace->element())),
  dofmap_ (functionSpace->dofmap()->copy()),
  subElements_ (2,nullptr),
  subDofmaps_ (2,nullptr)
  //dofmap_ (new dolfin::DofMap(std::shared_ptr<const ufc::dofmap>(new const poisson_dofmap_1), *mesh_))
  //dofmap_ (new dolfin::DofMap(*functionSpace->dofmap(),*mesh_))
  //dofmap_ (new const dolfin::DofMap(static_cast<const dolfin::DofMap&>(*functionSpace->dofmap())))
  {
    //dolfin::DofMapBuilder::build(*dofmap_,*mesh_,nullptr);
//    dofmap_ = dofmap_->create(*mesh_);
  }

template <class T_MeshMover, class T_FunctionSpace>
MeshManager<T_MeshMover,T_FunctionSpace>::MeshManager (const T_MeshMover& meshMover, T_FunctionSpace& functionSpace) :
  meshMover_ (meshMover),
  functionSpace_ (&functionSpace),
  mesh_ (new dolfin::Mesh(*functionSpace.mesh())),
  element_ (new const dolfin::FiniteElement(*functionSpace->element())),
  dofmap_ (functionSpace->dofmap()->copy()),
  subElements_ (2,nullptr),
  subDofmaps_ (2,nullptr)
  //dofmap_ (new const dolfin::DofMap(static_cast<const dolfin::DofMap&>(*functionSpace->dofmap())))
  {
    //dolfin::DofMapBuilder::build(*dofmap_,*mesh_,nullptr);
//    dofmap_ = dofmap_->create(*mesh_);
  }

template <class T_MeshMover, class T_FunctionSpace>
const T_FunctionSpace& MeshManager<T_MeshMover,T_FunctionSpace>::functionSpace()
  {
    return *functionSpace_;
  }

template<class T_MeshMover, class T_FunctionSpace>
void MeshManager<T_MeshMover,T_FunctionSpace>::moveMesh (const dolfin::GenericFunction& displacement)
{
  std::cout << "muovo la mesh di " << functionSpace_->str(false) << std::endl;
  meshMover_.move (*mesh_, displacement);
  *functionSpace_ = T_FunctionSpace(mesh_,element_,dofmap_);
//  for (auto iter=functionSpace_->component().begin(); iter != functionSpace_->component().end(); iter++)
//  {
//    auto comp((*functionSpace_)[*iter]);
  for (int i=0; i!=2; i++)
  {
//    std::cerr << "componente " << *iter << "  " << comp->str(false) << std::endl;
    std::cerr << "componente " << i << "  " << (*functionSpace_)[i]->str(false) << std::endl;
    std::string title("mesh della componente "); title += std::to_string(i);
    dolfin::plot(*(*functionSpace_)[i]->mesh(),title.c_str());
    dolfin::interactive();
    // sistema le componenti
    subElements_[i] = std::shared_ptr<dolfin::FiniteElement> (new dolfin::FiniteElement(*(*functionSpace_)[i]->element()));
    subDofmaps_[i]  = std::shared_ptr<dolfin::GenericDofMap>((*functionSpace_)[i]->dofmap()->copy());
    (*(*functionSpace_)[i]) = T_FunctionSpace(mesh_,subElements_[i],subDofmaps_[i]);
    title += " dopo";
    dolfin::plot(*(*functionSpace_)[i]->mesh(),title.c_str());
    dolfin::interactive();
  }

}

} // end of namespace

#endif
