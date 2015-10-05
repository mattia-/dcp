#ifndef MESHMANAGER_H_INCLUDE_GUARD
#define MESHMANAGER_H_INCLUDE_GUARD

#include <dcp/differential_problems/AbstractProblem.h>
#include <dcp/differential_problems/LinearProblem.h>
#include "laplaceVec.h"
#include "laplace1D.h"
#include "utilities.h"
#include <limits>

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
 		 *										1 : the boundary on which the given displacement has to be imposed
		 * 										2 : fixed boundary
		 *										3 : slip boundary
		 *										0 : interior facets / never mind
		 *	\param noSlipComponent the function component along which slip is not allowed,
		 *													on the boundary marked as 3 by meshFacets
		 *												 ( default value = std::numeric_limits<std::size_t>::max() )
		 */
  	MeshManager (const std::shared_ptr<dolfin::Mesh> mesh,
		 						 const std::shared_ptr<const dolfin::FacetFunction<std::size_t> > meshFacets,
	 						   const std::size_t noSlipComponent=std::numeric_limits<std::size_t>::max());
/*	  MeshManager (const std::shared_ptr<T_MeshMover> meshMover, std::shared_ptr<dolfin::FunctionSpace> functionSpace);
    MeshManager (const T_MeshMover& meshMover, dolfin::FunctionSpace& functionSpace);*/
    // !!! servono altri costruttori?

    void getDofOrder (const dolfin::FunctionSpace & funSp, std::size_t component=0);

//    const dolfin::FunctionSpace& functionSpace() const;
    std::shared_ptr<const dolfin::Mesh> mesh() const;
    
    void moveMesh ()
    {
        moveMesh (this->w_);
    }
    void moveMesh (const dolfin::Function& displacement, const std::string& component="all", const double dt=1.0);
//move//    std::shared_ptr<dolfin::MeshDisplacement> moveMesh (const dolfin::GenericFunction& boundaryDisplacement);
    std::shared_ptr <const dolfin::Function> computeDisplacement (const dolfin::Function& displacement, const std::string& component="all", const double dt=1.0);

    std::shared_ptr <const dolfin::Function> displacement () const;
    
    void updateMesh (std::shared_ptr<dolfin::FunctionSpace> functionSpace) const;

    //void setImposedDisplacement ();//const dolfin::GenericFunction & imposedDisplacement);
    void setImposedDisplacement (std::shared_ptr<const dolfin::FunctionSpace> imposedDisplacementSpace);

  protected:

    //void setMeshAll (std::shared_ptr<dolfin::FunctionSpace> functionSpace, std::shared_ptr<dolfin::Mesh> mesh);
    template <class T_FunctionSpace>
		    void setMeshAll (std::shared_ptr<T_FunctionSpace> functionSpace, std::shared_ptr<dolfin::Mesh> mesh);

    T_MeshMover meshMover_;
    const std::shared_ptr<dolfin::Mesh> mesh_;
    dolfin::FunctionSpace wFunSp_;
    dolfin::Function w_;
    bool displacementIsComputed_;
    std::shared_ptr<dcp::AbstractProblem> problemALE_;
/*    std::vector<std::shared_ptr<dolfin::FiniteElement> > subElements_;
    std::vector<std::shared_ptr<dolfin::GenericDofMap> > subDofmaps_;
*/
  	const std::shared_ptr<const dolfin::FacetFunction<std::size_t> > meshFacets_;
  	const std::size_t noSlipComponent_;
    dolfin::MeshFunction<bool> on_boundary_;
    dolfin::MeshFunction<bool> on_to_be_neglected_boundary_;
    std::vector<std::size_t> orderedMeshIdxs_;
    std::vector<std::size_t> orderedDisplacementDofs_;

    //ImposedDisplacement imposedDisplacement_;
    //ImposedFromFile imposedDisplacement_;
    ImposedDisplacement imposedDisplacement_;
    std::shared_ptr<const dolfin::FunctionSpace> imposedDisplacementSpace_;
    bool displacementIsImposed_;

};

// =============== IMPLEMENTATION ===============

template <class T_MeshMover>
MeshManager<T_MeshMover>::MeshManager (const std::shared_ptr<dolfin::Mesh> mesh,
 																			 const std::shared_ptr<const dolfin::FacetFunction<std::size_t> > meshFacets,
																			 const std::size_t noSlipComponent) :
    meshMover_ (),
		mesh_ (mesh),
		wFunSp_ (laplaceVec::FunctionSpace(*mesh_)),
		w_ (wFunSp_),
    displacementIsComputed_ (false),
		problemALE_ (new dcp::LinearProblem<laplaceVec::BilinearForm, laplaceVec::LinearForm> (dolfin::reference_to_no_delete_pointer(wFunSp_))),
		meshFacets_ (meshFacets),
		noSlipComponent_ (noSlipComponent),
		on_boundary_ (* mesh_, 0),
		on_to_be_neglected_boundary_ (* mesh_, 0),
    orderedMeshIdxs_ (mesh_->coordinates().size()/2),
    orderedDisplacementDofs_ (),
    displacementIsImposed_ (false)
{
    const std::vector<double> & coords (mesh_->coordinates());
//for (auto it=coords.begin(); it!=coords.end(); ++it) std::cerr << *it << ' '; std::cerr << std::endl;
    std::vector<std::tuple<std::size_t,double,double> > ordered;
    for (std::size_t k=0; k!=coords.size()/2; ++k)
      ordered.push_back (std::make_tuple (k, coords[2*k], coords[2*k+1]));
    std::sort (ordered.begin(),ordered.end(),compare);
    std::transform (ordered.begin(),ordered.end(),orderedMeshIdxs_.begin(),getFirst);
}

template <class T_MeshMover>
void MeshManager<T_MeshMover>::getDofOrder (const dolfin::FunctionSpace & funSp, std::size_t component)
{
    const std::vector<double> & dofsCoords (funSp.dofmap()->tabulate_all_coordinates (* mesh_));
    std::vector<dolfin::la_index> dofs (funSp[component]->dofmap()->dofs());
    std::vector<std::tuple<std::size_t,double,double,std::size_t> > ordered;
    for (std::size_t subcomp=0; subcomp!=funSp[component]->element()->value_dimension(0); ++subcomp)
    {
      std::vector<dolfin::la_index> dofs ((* funSp[component])[subcomp]->dofmap()->dofs());
      for (auto it=dofs.begin(); it!=dofs.end(); ++it)
      {
        ordered.push_back (std::make_tuple (*it, dofsCoords[2*(*it)], dofsCoords[2*(*it)+1], subcomp));
      }
    }
    std::sort (ordered.begin(),ordered.end(),compareDofs);
    orderedDisplacementDofs_.resize (ordered.size());
    std::transform (ordered.begin(),ordered.end(),orderedDisplacementDofs_.begin(),getFirstDofs);
//std::cerr << "orderedMeshIdxs_ (" << orderedMeshIdxs_.size() << ")   "; for (auto it=orderedMeshIdxs_.begin(); it!=orderedMeshIdxs_.end(); ++it) std::cerr << *it << ' '; std::cerr << std::endl;
//std::cerr << "orderedDisplacementDofs_ (" << orderedDisplacementDofs_.size() << ")   "; for (auto it=orderedDisplacementDofs_.begin(); it!=orderedDisplacementDofs_.end(); ++it) std::cerr << *it << ' '; std::cerr << std::endl;
//std::cerr << "dofsCoords (" << dofsCoords.size() << ")   "; for (auto it=dofsCoords.begin(); it!=dofsCoords.end(); ++it) std::cerr << *it << ' '; std::cerr << std::endl;
//std::cerr << " OK fino a " << __LINE__ << std::endl;
}

template <class T_MeshMover>
std::shared_ptr<const dolfin::Mesh> MeshManager<T_MeshMover>::mesh() const
{
    return mesh_;
}
													
template <class T_MeshMover>
void MeshManager<T_MeshMover>::moveMesh (const dolfin::Function& displacement, const std::string& component, const double dt)
{
    if (!displacementIsComputed_)
        computeDisplacement (displacement, component, dt);

 		// move the mesh and update mesh-dependent quantities
    mesh_->move (w_);
    //imposedDisplacement_.update ();
    displacementIsComputed_ = false;
}
    
template <class T_MeshMover>
std::shared_ptr <const dolfin::Function> MeshManager<T_MeshMover>::displacement () const
{
    return dolfin::reference_to_no_delete_pointer (w_);
}

template <class T_MeshMover>
void MeshManager<T_MeshMover>::updateMesh (std::shared_ptr<dolfin::FunctionSpace> functionSpace) const
{
 		 
}

template <class T_MeshMover>
template <class T_FunctionSpace>
void MeshManager<T_MeshMover>::setMeshAll (std::shared_ptr<T_FunctionSpace> functionSpace, std::shared_ptr<dolfin::Mesh> mesh)
{
// 		 TO BE DELETED
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
    displacementIsImposed_ = true;

    // avoid re-computation
    displacementIsComputed_ = true;

    return dolfin::reference_to_no_delete_pointer (w_);
  }

  const std::vector<double> & coords (mesh_->coordinates());
  mesh_->init (0,1);
//for (auto it=coords.begin(); it!=coords.end(); ++it) std::cerr << *it << ' '; std::cerr << std::endl;
//const std::vector<double> & funcoords (displacement[1].function_space()->mesh()->coordinates());
//for (auto it=funcoords.begin(); it!=funcoords.end(); ++it) std::cerr << *it << ' '; std::cerr << std::endl;
  std::vector<std::tuple<std::size_t,double,double> > ordered;
  for (std::size_t k=0; k!=coords.size()/2; ++k)
    ordered.push_back (std::make_tuple (orderedMeshIdxs_[k], coords[2*k], coords[2*k+1]));
  std::sort (ordered.begin(),ordered.end(),compareFirst);
//for (auto it=ordered.begin(); it!=ordered.end(); ++it) std::cerr << std::get<1>(*it) << ' ' << std::get<2>(*it) << '\t'; std::cerr << std::endl;
  std::vector<double> newCoords (coords.size());
  double x,y;
  std::size_t k, kDofs;
assert (ordered.size()==(problemData.nx+1)*(problemData.ny+1));
assert (orderedDisplacementDofs_.size()==2*((problemData.nx+1)*(problemData.ny+1)+3*problemData.nx*problemData.ny+problemData.nx+problemData.ny));
//assert (orderedDisplacementDofs_.size()==nx*ny+displacement.vector()->size());
  std::vector<double> newy (problemData.nx+1);
  for (std::size_t j=0; j!= newy.size(); ++j)
  {
    k = ordered.size()-newy.size()+j;
    kDofs = orderedDisplacementDofs_.size()-1 - 4*(newy.size()-1) + 4*j;
    //newy[j] = 0.1+std::get<1>(ordered[j])*(1-std::get<1>(ordered[j]));
    //newy[j] = 1+0.1*std::get<1>(ordered[j])+0.5*sin(2*3.1459*std::get<1>(ordered[j]));
  // vertical component of displacement
		/* newy[j] = std::get<2>(ordered[k])+dt*displacement[1](std::get<1>(ordered[k]),std::get<2>(ordered[k]));
    */
  // using dofs
    double val ( (* displacement.vector())[orderedDisplacementDofs_[kDofs]]);
		newy[j] = std::get<2>(ordered[k])+dt*val;
  // to ensure w.n = u.n
    /* x = std::get<1>(ordered[k]);
    y = std::get<2>(ordered[k])-1.97e-7;
//std::cerr << x << ' ' << y << ' '; std::cerr << displacement[1](x,y) << ' '; std::cerr << displacement[0](x,y) << '\t';
    if (j==0 || j==newy.size()-1)
		  newy[j] = std::get<2>(ordered[k])+dt*displacement[1](x,y);
    else
    {
      dolfin::Vertex v (* mesh_, k);
      std::array<double,2> vals, weights;
      std::size_t iii (0);
      for (dolfin::FacetIterator facet (v); !facet.end(); ++facet)
      {
        if (! facet->exterior())
          continue;
  	  	vals[iii] =  dt * ( facet->normal(0)*displacement[0](x,y) + facet->normal(1)*displacement[1](x,y) )
                      / facet->normal(1);
        weights[iii++] = facet->normal(1);
      }
      newy[j] = std::get<2>(ordered[k]) + ( vals[0] * weights[1] + vals[1]*weights[0] ) / (weights[0]+weights[1]);
    }
//std::cerr << std::get<1>(ordered[k]) << ' ' << std::get<2>(ordered[k]) << ' ' << displacement[1](std::get<1>(ordered[k]),std::get<2>(ordered[k])) << ' ' << newy[j] << std::endl;
    */
  }
//std::cerr << "displacement.vector   "; for (std::size_t iii=0; iii!=displacement.vector()->size(); ++iii) std::cerr << (*displacement.vector())[iii] << ' '; std::cerr << std::endl;
//std::cerr << "  last     " << std::setprecision(10) << newy[newy.size()-1] << std::endl;
//std::cerr << '\t' << k << '\t'; for (std::size_t iii=k; iii!=k-20; --iii) std::cerr << displacement[1](std::get<1>(ordered[iii]),std::get<2>(ordered[iii])) << ' '; std::cerr << std::endl;

//std::cerr << "newy    "; for (auto it=newy.begin(); it!=newy.end(); ++it) std::cerr << std::setprecision(10) << *it << ' '; std::cerr << std::endl;
assert (newy.back()!=0);
  for (std::size_t i=0; i!=problemData.ny+1; ++i)
  for (std::size_t j=0; j!=problemData.nx+1; ++j)
  {
    k = (problemData.nx+1)*i+j;
    x = std::get<1>(ordered[k]);
    y = std::get<2>(ordered[k]);
    newCoords [2*std::get<0>(ordered[k])] = x;
    //newCoords [2*std::get<0>(ordered[k])+1] = y + i/(double) ny * (newy[j]-std::get<2>(ordered[j]));
    newCoords [2*std::get<0>(ordered[k])+1] = i/(double) problemData.ny * newy[j];
  }
//std::cerr << "newCoords   "; for (auto it=newCoords.begin(); it!=newCoords.end(); ++it) std::cerr << *it << ' '; std::cerr << std::endl;
//std::cerr << '\t' << k << ' ' << x << ' ' << y << std::endl;

  dolfin::Mesh newmesh (* mesh_);
  newmesh.coordinates() = newCoords;
//dolfin::plot (newmesh,"newmesh"); dolfin::interactive();
  std::vector<double> newDofsCoords (wFunSp_.dofmap()->tabulate_all_coordinates (newmesh));
  std::vector<double> oldDofsCoords (wFunSp_.dofmap()->tabulate_all_coordinates (* mesh_));
  std::vector<double> dofDisplacement (oldDofsCoords.size());
  std::transform (newDofsCoords.begin(), newDofsCoords.end(), oldDofsCoords.begin(), dofDisplacement.begin(), std::minus<double>());
//std::cerr << "dofDisplacement   "; for (auto it=dofDisplacement.begin(); it!=dofDisplacement.end(); ++it) std::cerr << *it << ' '; std::cerr << std::endl;
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
//std::cerr << displacementVec.size() << ' ' << w_.vector()->size() << std::endl;
  displacementVec.apply("insert");
  * w_.vector() = displacementVec;

  // avoid re-computation
  displacementIsComputed_ = true;

  return dolfin::reference_to_no_delete_pointer (w_);
  
//  mesh.coordinates() = newCoords;
}
 								
template <class T_MeshMover>
//void MeshManager<T_MeshMover>::setImposedDisplacement ()//(const T_ImposedDisplacement & imposedDisplacement)
void MeshManager<T_MeshMover>::setImposedDisplacement (std::shared_ptr<const dolfin::FunctionSpace> imposedDisplacementSpace)
{
  displacementIsImposed_ = true;
  imposedDisplacementSpace_ = imposedDisplacementSpace;
}

#endif // MESHMANAGER_H_INCLUDE_GUARD

