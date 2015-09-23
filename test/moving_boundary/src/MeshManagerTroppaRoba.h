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
    std::shared_ptr <const dolfin::Function> computeDisplacementOld (const dolfin::Function& displacement, const std::string& component="all", const double dt=1.0);
    std::shared_ptr <const dolfin::Function> computeDisplacement (const dolfin::Function& displacement, const std::string& component="all", const double dt=1.0);

    std::shared_ptr <const dolfin::Function> displacement () const;
    
    void updateMesh (std::shared_ptr<dolfin::FunctionSpace> functionSpace) const;

    //void setImposedDisplacement ();//const dolfin::GenericFunction & imposedDisplacement);
    void setImposedDisplacement (std::shared_ptr<const dolfin::FunctionSpace> imposedDisplacementSpace);

  protected:

    //void setMeshAll (std::shared_ptr<dolfin::FunctionSpace> functionSpace, std::shared_ptr<dolfin::Mesh> mesh);
    template <class T_FunctionSpace>
		    void setMeshAll (std::shared_ptr<T_FunctionSpace> functionSpace, std::shared_ptr<dolfin::Mesh> mesh);
    void smooth (std::shared_ptr<dolfin::Mesh> mesh);
    dolfin::Function preSmooth (const dolfin::Function& displacement, const std::size_t component=1);

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
    dolfin::UnitIntervalMesh intervalMesh_;
    laplace1D::FunctionSpace intervalFunSp_;
    dolfin::Function intervalSol_;
    FunctionToExpression preSmoothed_;
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
		intervalMesh_ (ny),
		intervalFunSp_ (intervalMesh_),
		intervalSol_ (intervalFunSp_),
		preSmoothed_ (dolfin::reference_to_no_delete_pointer(intervalSol_),2,1),
    orderedMeshIdxs_ (mesh_->coordinates().size()/2),
    orderedDisplacementDofs_ (),
    displacementIsImposed_ (false)
{
    const std::vector<double> & coords (mesh_->coordinates());
for (auto it=coords.begin(); it!=coords.end(); ++it) std::cerr << *it << ' '; std::cerr << std::endl;
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
std::cerr << "orderedMeshIdxs_ (" << orderedMeshIdxs_.size() << ")   "; for (auto it=orderedMeshIdxs_.begin(); it!=orderedMeshIdxs_.end(); ++it) std::cerr << *it << ' '; std::cerr << std::endl;
std::cerr << "orderedDisplacementDofs_ (" << orderedDisplacementDofs_.size() << ")   "; for (auto it=orderedDisplacementDofs_.begin(); it!=orderedDisplacementDofs_.end(); ++it) std::cerr << *it << ' '; std::cerr << std::endl;
std::cerr << "dofsCoords (" << dofsCoords.size() << ")   "; for (auto it=dofsCoords.begin(); it!=dofsCoords.end(); ++it) std::cerr << *it << ' '; std::cerr << std::endl;
std::cerr << " OK fino a " << __LINE__ << std::endl;
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
std::shared_ptr <const dolfin::Function> MeshManager<T_MeshMover>::computeDisplacementOld (const dolfin::Function& displacement, const std::string& component, const double dt)
{
//std::cerr << wFunSp_.dofmap()->str(true) << std::endl;
//std::vector<double> coords (wFunSp_.dofmap()->tabulate_all_coordinates (* wFunSp_.mesh()));
//for (auto it=coords.begin(); it!=coords.end(); it++) std::cerr << *it << ' '; std::cerr << std::endl;
    // define SubDomain boundaries and basic coefficients
    dolfin::Constant zeroVec (0,0);
    dolfin::Constant zero (0);
    problemALE_->setCoefficient ("linear_form", dolfin::reference_to_no_delete_pointer (zeroVec), "zero");

 		// reset function spaces
 		
 		// set integration SubDomains in problemALE_
 		
 		// prepare the Dirichlet displacement datum
/*    dolfin::Function displ (* problemALE_->functionSpace());
    //imposedDisplacement_.setMesh (* mesh_);
    //imposedDisplacement_.loadFun ();
    //displ = imposedDisplacement_;
    FunctionToExpression preSmoothed (std::shared_ptr<dolfin::Function>(new dolfin::Function(preSmooth (displacement))),displacement.value_size(),1);
//    FunctionToExpressionBis preSmoothed (std::shared_ptr<dolfin::Function>(new dolfin::Function(preSmooth (displacement))),1,preSmoothed_.length());
    displ = preSmoothed;*/
 		
 		// remove + add DirichletBCs to problemALE_
std::unordered_map<std::size_t, double> bdval;
		problemALE_->removeDirichletBC ("imposedDisplacement");
 		//dolfin::DirichletBC dirBCimposed (* problemALE_->functionSpace(), displ, * this->meshFacets_, 1, "geometric");
 		dolfin::DirichletBC dirBCimposed (* (* problemALE_->functionSpace())[1], displacement[1], * this->meshFacets_, 1, "geometric");
dirBCimposed.get_boundary_values(bdval);
 		problemALE_->addDirichletBC (dirBCimposed, "imposedDisplacement");
		problemALE_->removeDirichletBC ("slipAllowed");
 		dolfin::DirichletBC dirBCslip (* (* problemALE_->functionSpace())[noSlipComponent_], zero, * this->meshFacets_, 3, "geometric");
dirBCslip.get_boundary_values(bdval);
 		problemALE_->addDirichletBC (dirBCslip, "slipAllowed");
/*		problemALE_->removeDirichletBC ("slipImposed");
 		dolfin::DirichletBC dirBCslip (* problemALE_->functionSpace(), displ, * this->meshFacets_, 3, "geometric");
dirBCslip.get_boundary_values(bdval);
 		problemALE_->addDirichletBC (dirBCslip, "slipImposed");*/
		problemALE_->removeDirichletBC ("fixedBoundary");
 		dolfin::DirichletBC dirBCfixed (* problemALE_->functionSpace(), zeroVec, * this->meshFacets_, 2, "geometric");
dirBCfixed.get_boundary_values(bdval);
 		problemALE_->addDirichletBC (dirBCfixed, "fixedBoundary");
 		
 		// solve problemALE_ and finalize solution
 		problemALE_->solve();
 		w_ = problemALE_->solution();
    (* w_.vector()) *= dt;
//this->preSmoothed_.updateRef(dt);
 		
//std::cerr << wFunSp_.dofmap()->str(true) << std::endl;
//coords = wFunSp_.dofmap()->tabulate_all_coordinates (* wFunSp_.mesh());
//for (auto it=coords.begin(); it!=coords.end(); it++) std::cerr << *it << ' '; std::cerr << std::endl;

    // avoid re-computation
    displacementIsComputed_ = true;

    return dolfin::reference_to_no_delete_pointer (w_);
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
void MeshManager<T_MeshMover>::smooth (std::shared_ptr<dolfin::Mesh> mesh)
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
dolfin::plot(*mesh, "post smoothing");dolfin::interactive();
    }
    mesh->smooth(5);
dolfin::plot(*mesh, "post dolfin smoothing");dolfin::interactive();
}

template <class T_MeshMover>
dolfin::Function MeshManager<T_MeshMover>::preSmooth (const dolfin::Function& displacement, const std::size_t component)
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
dolfin::plot(displ,* intervalSol_.function_space()->mesh(),"displ in preSmooth");dolfin::interactive();
    dolfin::DirichletBC bc2 (* intervalSol_.function_space(),displ,secondExtr);
    std::vector<const dolfin::DirichletBC*> bcs ({& bc2, & bc1});
        // l'ordine conta in questo vector!
        
    dolfin::solve (a==F,intervalSol_,bcs);
dolfin::plot(intervalSol_,"intervalSol");dolfin::interactive();

std::cerr << "salut" << intervalSol_.value_size() << std::endl;
//    FunctionToExpression returnSol (intervalSol,displacement.value_size(),component);

    return intervalSol_;//returnSol;
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
assert (ordered.size()==(nx+1)*(ny+1));
assert (orderedDisplacementDofs_.size()==2*((nx+1)*(ny+1)+3*nx*ny+nx+ny));
//assert (orderedDisplacementDofs_.size()==nx*ny+displacement.vector()->size());
  std::vector<double> newy (nx+1);
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
std::cerr << "  last     " << std::setprecision(10) << newy[newy.size()-1] << std::endl;
//std::cerr << '\t' << k << '\t'; for (std::size_t iii=k; iii!=k-20; --iii) std::cerr << displacement[1](std::get<1>(ordered[iii]),std::get<2>(ordered[iii])) << ' '; std::cerr << std::endl;

//std::cerr << "newy    "; for (auto it=newy.begin(); it!=newy.end(); ++it) std::cerr << std::setprecision(10) << *it << ' '; std::cerr << std::endl;
assert (newy.back()!=0);
  for (std::size_t i=0; i!=ny+1; ++i)
  for (std::size_t j=0; j!=nx+1; ++j)
  {
    k = (nx+1)*i+j;
    x = std::get<1>(ordered[k]);
    y = std::get<2>(ordered[k]);
    newCoords [2*std::get<0>(ordered[k])] = x;
    //newCoords [2*std::get<0>(ordered[k])+1] = y + i/(double) ny * (newy[j]-std::get<2>(ordered[j]));
    newCoords [2*std::get<0>(ordered[k])+1] = i/(double) ny * newy[j];
  }
//std::cerr << "newCoords   "; for (auto it=newCoords.begin(); it!=newCoords.end(); ++it) std::cerr << *it << ' '; std::cerr << std::endl;
std::cerr << '\t' << k << ' ' << x << ' ' << y << std::endl;

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
  for (std::size_t k=0; k!=processedDofDisplacement.size()/2; ++k)
  {
    processedDofDisplacement [2*k] = dofDisplacement[4*k];
    processedDofDisplacement [2*k+1] = dofDisplacement[4*k+1];
  }
  dolfin::Vector displacementVec (MPI_COMM_WORLD, processedDofDisplacement.size());
  displacementVec.set_local (processedDofDisplacement);
std::cerr << displacementVec.size() << ' ' << w_.vector()->size() << std::endl;
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

