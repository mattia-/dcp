#include <dolfin.h>
#include "laplaceVec.h"
#include <string>

class MyBoundary : public dolfin::SubDomain
{
  public :
    bool inside (const dolfin::Array<double>& x, bool on_boundary) const
    {
      return dolfin::near(x[0],0) || dolfin::near(x[1],0);
    }
};

class Riletto : public dolfin::Expression
{
  public :
    void setFunSp (dolfin::FunctionSpace& funSp)
    {
        funSp_.reset (&funSp);
    }

    void loadFun (std::string funname, std::string filename)
    {
        file_.reset (new dolfin::HDF5File(MPI_COMM_WORLD,filename,"r"));
        fun_.reset (new dolfin::Function(*funSp_));
        file_->read(*fun_,funname);
    }
  
    void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x) const
    {
        (* fun_)[1].eval (values,x);
    }
/*        fun_->eval (values,x);
    }
    std::size_t value_rank() const
    {
      return 1;
    }

    std::size_t value_dimension(std::size_t i) const
    {
      return 2;
    }
*/

  private :
    std::shared_ptr<dolfin::HDF5File> file_;
    std::shared_ptr<dolfin::FunctionSpace> funSp_;
    std::shared_ptr<dolfin::Function> fun_;
    
};

int main()
{
  dolfin::RectangleMesh mesh (0,0,1,2,16,32,"left");
  laplaceVec::FunctionSpace fs (mesh);
  
  dolfin::FacetFunction<std::size_t> meshFun (mesh);
  meshFun.set_all(0);
  MyBoundary bd;
  bd.mark(meshFun,1);

  dolfin::plot(meshFun,"meshFun"); dolfin::interactive();
  dolfin::Constant values (0,3);
  std::unordered_map<std::size_t, double> valMap;

  dolfin::DirichletBC dirBCtop (fs,values,meshFun,1,"topological");
  dolfin::DirichletBC dirBCgeo (fs,values,meshFun,1,"geometric");
  dolfin::DirichletBC dirBCptw (fs,values,meshFun,1,"pointwise");
  dirBCtop.get_boundary_values(valMap); std::cerr << valMap.size() << std::endl;
  dirBCgeo.get_boundary_values(valMap); std::cerr << valMap.size() << std::endl;
//  dirBCptw.get_boundary_values(valMap); std::cerr << valMap.size() << std::endl;

  dolfin::DirichletBC dirBCtop2 (fs,values,bd,"topological");
  dolfin::DirichletBC dirBCgeo2 (fs,values,bd,"geometric");
  dolfin::DirichletBC dirBCptw2 (fs,values,bd,"pointwise");
  dirBCtop2.get_boundary_values(valMap); std::cerr << valMap.size() << std::endl;
  dirBCgeo2.get_boundary_values(valMap); std::cerr << valMap.size() << std::endl;
//  dirBCptw2.get_boundary_values(valMap); std::cerr << valMap.size() << std::endl;

  dolfin::Function u (fs);
  laplaceVec::BilinearForm a (fs,fs);
  laplaceVec::LinearForm L (fs);

  std::shared_ptr<dolfin::Constant> zero (new dolfin::Constant(0,0));
  L.set_coefficient("zero",zero);

  dolfin::solve (a==L, u, dirBCgeo);
  dolfin::plot (u); dolfin::interactive();

  dolfin::HDF5File uOutFile (MPI_COMM_WORLD,"prova.bin","w");
//  uFile << std::make_pair(&(u[1]),4e-5);
  uOutFile.write (u[1],"velY");
  uOutFile.write (u[0],"velX");
  uOutFile.write (u,"velAll");
  uOutFile.close();

  dolfin::HDF5File uInFile (MPI_COMM_WORLD,"prova.bin","r");
  dolfin::Function u1 (fs[1]->collapse());//,"prova.bin");
  uInFile.read(u1,"velY");
//  dolfin::plot (u1); dolfin::interactive();
  dolfin::Function u0 (fs[0]->collapse());//,"prova.bin");
  uInFile.read(u0,"velX");
//  dolfin::plot (u0); dolfin::interactive();
  dolfin::Function uR (fs);//,"prova.bin");
  uInFile.read(uR,"velAll");
//  dolfin::plot (uR); dolfin::interactive();
  uInFile.close();

  Riletto uRil;
  uRil.setFunSp (fs);
  uRil.loadFun ("velAll","prova.bin");
  dolfin::plot (uRil,mesh); dolfin::interactive();

  return 0;
}
