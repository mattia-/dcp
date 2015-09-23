#include <dolfin.h>
#include "utilities.h"
#include <vector>
// provaA
/*#include "laplace1D.h"*/
// provaB provaC
#include "laplaceVec.h"
#include "myNavierstokesTimeCurvLinear.h"
#include <functional>
#include <algorithm>

// provaA
/*dolfin::Function bla ()
{
  dolfin::UnitIntervalMesh mesh (5);
  laplace1D::FunctionSpace fs (mesh);
  dolfin::Function intervalSol (fs);
  FunctionToExpression blabla (intervalSol);
  return intervalSol;//blabla;
}*/

// provaB
/*class Displacement : public dolfin::Expression
{
    void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x) const
    {
        values[0] = 0.0;
        values[1] = x[1];
    }
    
    std::size_t value_rank() const
    {
      return 1;
    }

    std::size_t value_dimension(std::size_t i) const
    {
      return 2;
    }
};*/

// provaC
/*bool compare (const std::tuple<std::size_t,double,double> p1, const std::tuple<std::size_t,double,double> p2)
{
  return ( ( std::get<2>(p1) < std::get<2>(p2) ) || ( std::get<2>(p1) == std::get<2>(p2) && std::get<1>(p1) < std::get<1>(p2)) );
}*/

int main()
{
// provaA
/*  dolfin::UnitIntervalMesh mesh (5);
  laplace1D::FunctionSpace fs (mesh);
  dolfin::Function fun (fs);
  FunctionToExpression a (fun);
  FunctionToExpression b (a);
  std::cerr << b.value_size() << ' ' << a.value_size() << std::endl;
  std::cerr << b(0.2) << ' ' << a(0.2) << std::endl;

  FunctionToExpression c (bla());
  std::cerr << c.value_size() << c(0.2) << std::endl;*/

// provaB
/*dolfin::parameters["allow_extrapolation"] = true;
  dolfin::UnitSquareMesh mesh (10,10);
  laplaceVec::FunctionSpace fs (dolfin::reference_to_no_delete_pointer(mesh));
  dolfin::Function fun (fs);
  fun = (XY());
  dolfin::plot (fun,"initial"); dolfin::interactive();
  mesh.move (Displacement());
  dolfin::plot (fun,"moved"); dolfin::interactive();
  fun = (XY());
  dolfin::plot (fun,"re-computed"); dolfin::interactive();*/

// provaC
  std::size_t nx (2), ny (3);
  dolfin::UnitSquareMesh mesh (nx,ny);
  dolfin::plot(mesh, "initial mesh"); dolfin::interactive();
  const std::vector<double> & coords (mesh.coordinates());
for (auto it=coords.begin(); it!=coords.end(); ++it) std::cerr << *it << ' '; std::cerr << std::endl;
  std::vector<std::tuple<std::size_t,double,double> > ordered;
  for (std::size_t k=0; k!=coords.size()/2; ++k)
    ordered.push_back (std::make_tuple (k, coords[2*k], coords[2*k+1]));
  std::sort (ordered.begin(),ordered.end(),compare);
for (auto it=ordered.begin(); it!=ordered.end(); ++it) std::cerr << std::get<1>(*it) << ' ' << std::get<2>(*it) << '\t'; std::cerr << std::endl;
  std::vector<double> newCoords (coords.size());
  double x,y;
  std::size_t k;
assert (ordered.size()==(nx+1)*(ny+1));
  std::vector<double> newy (nx+1);
  for (std::size_t j=0; j!= newy.size(); ++j)
    //newy[j] = 0.1+std::get<1>(ordered[j])*(1-std::get<1>(ordered[j]));
    newy[j] = 1+0.1*std::get<1>(ordered[j])+0.5*sin(2*3.1459*std::get<1>(ordered[j]));
  for (std::size_t i=0; i!=ny+1; ++i)
  for (std::size_t j=0; j!=nx+1; ++j)
  {
    k = (nx+1)*i+j;
    x = std::get<1>(ordered[k]);
    y = std::get<2>(ordered[k]);
std::cerr << x << ' ' << y << std::endl;
    newCoords [2*std::get<0>(ordered[k])] = x;
    //newCoords [2*std::get<0>(ordered[k])+1] = y;
    //newCoords [2*std::get<0>(ordered[k])+1] = std::get<2>(ordered[j]) + i/(double) ny * (std::get<2>(ordered[(nx+1)*ny+j])-std::get<2>(ordered[j]));
    newCoords [2*std::get<0>(ordered[k])+1] = std::get<2>(ordered[j]) + i/(double) ny * (newy[j]-std::get<2>(ordered[j]));
  }
std::cerr << "newCoords   "; for (auto it=newCoords.begin(); it!=newCoords.end(); ++it) std::cerr << *it << ' '; std::cerr << std::endl;
  dolfin::Mesh newmesh (mesh);
  dolfin::Mesh oldmesh (mesh);
  newmesh.coordinates() = newCoords;
  dolfin::plot(newmesh, "new mesh"); dolfin::interactive();
  laplaceVec::FunctionSpace fs (mesh);
  myNavierstokesTimeCurvLinear::FunctionSpace fsNS (mesh);
std::cerr << "spaces:                fs dimension   " << fs.dim() << '\t' << "fsNS dimension   " << fsNS.dim() << '\t' << "fsNS[0] dimension   " << fsNS[0]->dim() << std::endl ;
std::cerr << "elements(space_dim):   fs dimension   " << fs.element()->space_dimension() << '\t' << "fsNS dimension   " << fsNS.element()->space_dimension() << '\t' << "fsNS[0] dimension   " << fsNS[0]->element()->space_dimension() << std::endl ;
std::cerr << "elements(value_dim):   fs dimension   " << fs.element()->value_dimension(0) << '\t' << "fsNS dimension   " << fsNS.element()->value_dimension(0) << '\t' << "fsNS[0] dimension   " << fsNS[0]->element()->value_dimension(0) << std::endl ;
//  laplaceVec::FunctionSpace oldfs (oldmesh);
  dolfin::Function w (fs);
//  dolfin::Function oldw (fs);
  const dolfin::GenericDofMap& dofmap (* fs.dofmap());
  const dolfin::GenericDofMap& dofmapNS (* fsNS.dofmap());
  const dolfin::GenericDofMap& dofmapNSVel (* fsNS[0]->dofmap());
  const dolfin::GenericDofMap& dofmapNSVel0 (* (* fsNS[0])[0]->dofmap());
  const dolfin::GenericDofMap& dofmapNSVel1 (* (* fsNS[0])[1]->dofmap());
  //!!! la riga seguente induce segfault quando si fa dofmapNSVel.dofs()
  //const dolfin::GenericDofMap& dofmapNSVel (* dofmapNS.extract_sub_dofmap ({0}, * fsNS.mesh()));
  std::vector<dolfin::la_index> dofs (dofmap.dofs());
  std::vector<dolfin::la_index> dofsNS (dofmapNS.dofs());
  std::vector<dolfin::la_index> dofsNSVel (dofmapNSVel.dofs());
  std::vector<dolfin::la_index> dofsNSVel0 (dofmapNSVel0.dofs());
  std::vector<dolfin::la_index> dofsNSVel1 (dofmapNSVel1.dofs());
std::cerr << "velocity dofs   "; for (auto it=dofsNSVel.begin(); it!=dofsNSVel.end(); ++it) std::cerr << *it << ' '; std::cerr << std::endl;
std::cerr << "x velocity dofs   "; for (auto it=dofsNSVel0.begin(); it!=dofsNSVel0.end(); ++it) std::cerr << *it << ' '; std::cerr << std::endl;
std::cerr << "y velocity dofs   "; for (auto it=dofsNSVel1.begin(); it!=dofsNSVel1.end(); ++it) std::cerr << *it << ' '; std::cerr << std::endl;
  std::vector<double> newDofsCoords (dofmap.tabulate_all_coordinates (newmesh));
std::cerr << "newDofsCoords   "; for (auto it=newDofsCoords.begin(); it!=newDofsCoords.end(); ++it) std::cerr << *it << ' '; std::cerr << std::endl;
//  const dolfin::GenericDofMap& oldDofmap (* oldfs.dofmap());
//  std::vector<dolfin::la_index> oldDofs (oldDofmap.dofs());
//  std::vector<double> oldDofsCoords (oldDofmap.tabulate_all_coordinates (oldmesh));
  std::vector<double> oldDofsCoords (dofmap.tabulate_all_coordinates (oldmesh));
std::cerr << "oldDofsCoords   "; for (auto it=oldDofsCoords.begin(); it!=oldDofsCoords.end(); ++it) std::cerr << *it << ' '; std::cerr << std::endl;
  std::vector<double> oldDofsCoordsNS (dofmapNS.tabulate_all_coordinates (oldmesh));
std::cerr << "oldDofsCoordsNS   "; for (auto it=oldDofsCoordsNS.begin(); it!=oldDofsCoordsNS.end(); ++it) std::cerr << *it << ' '; std::cerr << std::endl;
  std::vector<double> dofDisplacement (oldDofsCoords.size());
  std::transform (newDofsCoords.begin(), newDofsCoords.end(), oldDofsCoords.begin(), dofDisplacement.begin(), std::minus<double>());
std::cerr << "dofDisplacement   "; for (auto it=dofDisplacement.begin(); it!=dofDisplacement.end(); ++it) std::cerr << *it << ' '; std::cerr << std::endl;
  std::vector<double> processedDofDisplacement (dofDisplacement.size()/2);
    // NB assumiamo che in dofmap.dofs() i dof siano sempre messi a coppie, con la componente x prima e quella y dopo (quindi le coordinate si ripetono)
  for (std::size_t k=0; k!=processedDofDisplacement.size()/2; ++k)
  {
    processedDofDisplacement [2*k] = dofDisplacement[4*k];
    processedDofDisplacement [2*k+1] = dofDisplacement[4*k+1];
  }
  dolfin::Vector displacementVec (MPI_COMM_WORLD, processedDofDisplacement.size());
  displacementVec.set_local (processedDofDisplacement);
std::cerr << displacementVec.size() << ' ' << w.vector()->size() << std::endl;
  displacementVec.apply("insert");
  * w.vector() = displacementVec;
  dolfin::plot(w, "displacement on old mesh"); dolfin::interactive();
  mesh.coordinates() = newCoords;
  dolfin::plot(w, "displacement on new mesh"); dolfin::interactive();

  return 0;
}
