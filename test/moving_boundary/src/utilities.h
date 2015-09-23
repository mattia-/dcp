/* 
 *  Copyright (C) 2015, Ivan Fumagalli, ivan.fumagalli@polimi.it
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

#ifndef IVAN_UTILITIES_H_INCLUDE_GUARD
#define IVAN_UTILITIES_H_INCLUDE_GUARD

//#define GerbeauLelievre
//#define SprittlesShikhmurzaev
#define Yamamoto

//#define vanMourik

#include <dolfin.h>
#include <string>
#include <math.h>
#include "myNavierstokesTimeCurvLinear.h"

// parameters (defined in utilities.cpp)
extern std::string savepath;
extern double lx,ly,yxratio,ustar,TP_EPS;
extern std::size_t nx,ny;
extern double re,st,ca;
extern double betaVal,cosThetaSVal;
#ifdef Yamamoto
extern double rho;
#endif

class MovingLid : public dolfin::SubDomain
{
    bool inside (const dolfin::Array<double>& x, bool on_boundary) const
    {
        return dolfin::near (x[1], 1) && !(dolfin::near (x[0],0)) && on_boundary;
    }
};

class FixedWalls : public dolfin::SubDomain
{
    bool inside (const dolfin::Array<double>& x, bool on_boundary) const
    {
        return (dolfin::near (x[1], 0) && on_boundary)
               ||
               (dolfin::near (x[0], 0) && on_boundary) 
               ||
               (dolfin::near (x[0], 1) && on_boundary);
               
    }
};
class DirichletBoundary : public dolfin::SubDomain
{
    bool inside (const dolfin::Array<double>& x, bool on_boundary) const
    {
        return  on_boundary
                && (
                dolfin::near(x[1], 0)
                ||
                dolfin::near(x[0], 0)
                ||
                dolfin::near(x[0], 1) );
    }
};
class DirichletData : public dolfin::Expression
{
public:

  void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x) const
  {
    values[0]= 0.0;// * x[0] * (1-x[0]);
    values[1]= 0.2 * x[0] * (1-x[0]) * (0.5*lx-x[0]);
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
class DirichletDataComplete : public dolfin::Expression
{
public:

  void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x) const
  {
    values[0]= 0.0;// * x[0] * (1-x[0]);
    values[1]= 0.01 * x[0] * (1-x[0]) * x[1];
    values[2]= 0.0;
  }

  std::size_t value_rank() const
  {
    return 1;
  }

  std::size_t value_dimension(std::size_t i) const
  {
    return 3;
  }

};
class LateralInterior : public dolfin::SubDomain
{
    bool inside (const dolfin::Array<double>& x, bool on_boundary) const
    {
        return on_boundary
               && ( dolfin::near (x[0],0.0) || dolfin::near(x[0],lx) )
               && ( x[1] < ly-6e-16);
    }

    friend class TriplePoints;
    friend class TriplePointLeft;
    friend class TriplePointRight;
};
class LateralWall : public dolfin::SubDomain
{
  public:
    bool inside (const dolfin::Array<double>& x, bool on_boundary) const
    {
        return on_boundary
               && ( dolfin::near (x[0],0.0,0.5*TP_EPS) || dolfin::near(x[0],lx) );
    }
};
class TopBoundary : public dolfin::SubDomain
{
  public :
    bool inside (const dolfin::Array<double>& x, bool on_boundary) const
    {
        return on_boundary
               &&
//              ( x[1] >= ly - 6.0e-16 );
               dolfin::near (x[1],ly);
/*          //what follows takes away the triple points
               &&
             !( dolfin::near(x[0],0) )
               &&
             !( dolfin::near(x[0],lx) );*/
    }

  private:
    friend class TriplePointLeftVertex;
    friend class TriplePointRightVertex;
};
class BottomBoundary : public dolfin::SubDomain
{
  public:
    bool inside (const dolfin::Array<double>& x, bool on_boundary) const
    {
        return on_boundary
               &&
              (dolfin::near (x[1],0.0));
    }
};
/*class TriplePoints : public dolfin::SubDomain
{
    bool inside (const dolfin::Array<double>& x, bool on_boundary) const
    {
        return on_boundary
              && ( ( dolfin::near(x[0],0,0.006) && dolfin::near(x[1],ly,0.006) )
                ||
               ( dolfin::near(x[0],lx,0.006) && dolfin::near(x[1],ly,0.006) ) )
              && !( lateralInterior.inside(x, on_boundary) );
        return  ( dolfin::near(x[0], 0) && (x[1]<0.01-6e-16) && (x[1]>6e-16)) //dolfin::near(x[1],0.005) )
       		     ||
             		( dolfin::near(x[0], 0.05) && (x[1]<0.01-6e-16) && (x[1]>6e-16)); //dolfin::near(x[1],0.005) );
    }

    LateralInterior lateralInterior;
};*/
class TriplePointLeft : public dolfin::SubDomain
{
    bool inside (const dolfin::Array<double>& x, bool on_boundary) const
    {
        return on_boundary
              && ( dolfin::near(x[0],0,TP_EPS) && dolfin::near(x[1],ly) )
              && !( lateralInterior.inside(x, on_boundary) );
    }

    LateralInterior lateralInterior;
};
class TriplePointRight : public dolfin::SubDomain
{
    bool inside (const dolfin::Array<double>& x, bool on_boundary) const
    {
        return on_boundary
              && ( dolfin::near(x[0],lx,TP_EPS) && dolfin::near(x[1],ly) )
              && !( lateralInterior.inside(x, on_boundary) );
    }

    LateralInterior lateralInterior;
};
class TriplePointLeftVertex : public dolfin::SubDomain
{
  public:
    bool inside (const dolfin::Array<double>& x, bool on_boundary) const
    {
        return on_boundary
              //&& ( dolfin::near(x[0],0) && top_.inside(x,on_boundary) );
              && ( dolfin::near(x[0],0) && dolfin::near(x[1],ly) );
    }

  private:
    TopBoundary top_;
};
class TriplePointRightVertex : public dolfin::SubDomain
{
  public:
    bool inside (const dolfin::Array<double>& x, bool on_boundary) const
    {
        return on_boundary
              //&& ( dolfin::near(x[0],lx) && top_.inside(x,on_boundary) );
              && ( dolfin::near(x[0],lx) && dolfin::near(x[1],ly) );
    }

  private:
    TopBoundary top_;
};

//TODO
/*class TriplePointVertices : public dolfin::SubDomain
{
  // To be used with the option "pointwise" in dolfin::DirichletBC
  // NB: "pointwise" gives always False as value for on_boundary 
  public:
    bool inside (const dolfin::Array<double>& x, bool on_boundary) const
    {
        return SERVE_UNA_MESHFUNCTION...
    }

  private:
  TriplePointLeftVertex leftVert_;
  TriplePointRightVertex rightVert_;
};*/

class InflowDirichletData : public dolfin::Expression
{
    void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x) const
    {
        values[0] = 0.0;
        values[1] = ustar * 4/(lx*lx)*x[0]*(lx-x[0]) ;
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

class XY : public dolfin::Expression
{
    void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x) const
    {
        values[0] = x[0];
        values[1] = x[1] ;
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
class XYZ : public dolfin::Expression
{
    void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x) const
    {
        values[0] = x[0];
        values[1] = x[1] ;
        values[2] = 0;
    }
    
    std::size_t value_rank() const
    {
      return 1;
    }

    std::size_t value_dimension(std::size_t i) const
    {
      return 3;
    }
};

class BottomBoundaryEvaluator
{ 
    public:
        bool operator() (const dolfin::Array<double>& x, bool on_boundary)
        {
            return dolfin::near (x[1], 0) && on_boundary;
        }
}; 
class InflowDirichletBCEvaluator
{
    public:
        void operator() (dolfin::Array<double>& values, const dolfin::Array<double>& x, const double& t)
        {
            values[0] = 0;
//            values[1] = sin(2*3.14*t) * 6*x[0]*(1-x[0]); //firstDrop
//            values[1] = 0.5 * sin(2*3.14*t) * 4*x[0]*(1-x[0]); //freq1
            values[1] = ustar * sin(t* 2*3.14) * 4/(lx*lx)*x[0]*(lx-x[0]);
        }
};

class BoundaryStressEvaluator
{

    public:
        BoundaryStressEvaluator () = delete;
        BoundaryStressEvaluator (double nominalStressFirst, double nominalStressSecond) :
          nominalStress ({nominalStressFirst,nominalStressSecond})
          {}

        void operator() (dolfin::Array<double>& values, const dolfin::Array<double>& x, const double& t)
        {
            values[0] = 0;
            values[1] = 9.81 * (ly-x[1]) * nominalStress.second;
        }

    private:
        std::pair<double,double> nominalStress;
};

class InitialDisplacement : public dolfin::Expression
{
public:

  void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x) const
  {
//      values[0] = -x[0] + ( x[0]<0.5*lx ? pow(x[0],a)*pow(0.5*lx,1-a) : (x[0]>0.5*lx ? lx-pow(lx-x[0],a)*pow(0.5*lx,1-a) : x[0]) );
//      values[1] = -x[1] + ly-pow(ly-x[1],b)*pow(ly,1-b);
      values[0] = 0;//-x[0] + ( x[0]<0.5*lx ? (exp(a*x[0])-1)*0.5*lx/(exp(a*0.5*lx)-1) : (x[0]>0.5*lx ? lx-(exp(a*(lx-x[0]))-1)*0.5*lx/(exp(a*0.5*lx)-1) : x[0]) );
      values[1] = -x[1] + ly-(exp(b*(ly-x[1]))-1)*ly/(exp(b*ly)-1);
  }

  std::size_t value_rank() const
  {
    return 1;
  }

  std::size_t value_dimension(std::size_t i) const
  {
    return 2;
  }

private:
  double a=5;
  double b=2e3;

};

class WallVelocity : public dolfin::Expression
{
public:

  void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x) const
  {
      values[0] = 0;
#if defined(GerbeauLelievre) && !defined(SprittlesShikhmurzaev) && !defined(Yamamoto)
      values[1] = x[0]<0.5*lx ? -0.25 : 0.25;
#elif defined(SprittlesShikhmurzaev) && !defined(GerbeauLelievre) && !defined(Yamamoto)
      values[1] = -1;
#elif defined(Yamamoto) && !defined(GerbeauLelievre) && !defined(SprittlesShikhmurzaev)
      values[1] = 0;
#endif
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

/*class MeshFunctionWrapper
{
  public:
    MeshFunctionWrapper(dolfin::FacetFunction<std::size_t> meshFacets) :
        meshFacets_(meshFacets)
    {}

    bool inside (const dolfin::Array<double>& x, bool on_boundary) const
    {
        return ...
    }

  private:
    dolfin::FacetFunction<std::size_t> meshFacets_;
};

class MeshFunctionToSubDomain : public dolfin::SubDomain
{
public :

    bool inside (const dolfin::Array<double>& x, bool on_boundary) const
    {
        return meshFunWrapper_.inside(x,on_boundary);
    }

    MeshFunctionWrapper meshFunWrapper_;
};*/

class ALEStiffnessVec : public dolfin::Expression
{
public:

  void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x) const
  {
      values[0] = 1;
#if defined(Yamamoto) && !defined(GerbeauLelievre) && !defined(SprittlesShikhmurzaev)
      values[1] = 1e3;//1+40*exp(x[1]/ly);
#endif
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

class myDouble
{
  // serve così sono sicuro di inizializzare a quello che voglio senza dover ridefinire il costruttore della classe che contiene un membro di tipo myDouble
    public:
//    myDouble() : initVal(0.000001), val(0.000001) {}
//    myDouble() : initVal(0.05), val(0.05) {}
    myDouble() : initVal(1), val(1) {}
    double initVal;
    mutable double val;
};
class myInt
{
  public:
  myInt() : val(0) {};
  
  std::size_t val;
};
/*class ImposedDisplacement : public dolfin::Expression
{
public:

  void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x) const
  {
      values[0] = 0;
//      values[1] = a.val * (x[0]-0.5*lx)*(x[0]-0.5*lx)/(0.5*lx*0.5*lx)/dt;
      values[1] = ( a.val>0 ) ? 
//                              (a.val/dt * 0.5*(1.1+cos(2.0*3.14159265/lx*x[0])))
//                              (a.val/dt * (0.05 + (x[0]-0.5*lx)*(x[0]-0.5*lx)/(0.5*lx*0.5*lx)))
                              (0.1*a.initVal + a.val * ( 0.05 + (exp(pow((x[0]-0.5*lx)/(0.5*lx),8))-1) / (exp(1)-1) ))
                              :
                              (0.1*a.initVal - 0.4*a.val * (0.05 + 1 - (x[0]-0.5*lx)*(x[0]-0.5*lx)/(0.5*lx*0.5*lx)));
//                              (0.1*a.initVal - 0.4*a.val * ( 0.05 + 1 - (exp(pow((x[0]-0.5*lx)/(0.5*lx),2))-1) / (exp(1)-1) ));
      values[1] = ( a.val > 0 ) ?
                                (bar + a.val * pow(x[0]-0.5*lx,4)/pow(0.5*lx,4) + (a.initVal-a.val) * (x[0]-lx/3.0)*(x[0]-lx/3.0)/(lx*lx/9.0)*(x[0]-2.0*lx/3.0)*(x[0]-2.0*lx/3.0)/(4.0*lx*lx/9.0))
                                :
                                (bar + b.val * (x[0]-lx/3.0)*(x[0]-lx/3.0)/(lx*lx/9.0)*(x[0]-2.0*lx/3.0)*(x[0]-2.0*lx/3.0)/(4.0*lx*lx/9.0) + (b.initVal-b.val) * (1-(x[0]-0.5*lx)*(x[0]-0.5*lx)/(lx*lx/4.0)));

      double time = timeStep.val*dt;
      double timeWeight = 9.79e+08*pow(time,4) - 1.08e+07*pow(time,3) + 4.26e+04*time*time - 70.9*time + 0.0700;
      values[0] *= timeWeight;
      values[1] *= timeWeight;
  }

  std::size_t value_rank() const
  {
    return 1;
  }

  std::size_t value_dimension(std::size_t i) const
  {
    return 2;
  }

  void update()
  {
//      a.val -= a.initVal/nT;
//      b.val -= (a.val>0) ? 0 : b.initVal/(2*nT);
      ++(timeStep.val);
  }
private:
  myDouble a,b;
  myInt timeStep;
//  double b = 0.0000002;
  double bar = 0.03;
  std::size_t nT = 100;
  double dt = 0.00002;

};*/
class ImposedDisplacement : public dolfin::Expression
{
public:

  void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x) const
  {
      values[0] = 0;
      values[1] = 0.05;
      values[2] = 0;
  }

  std::size_t value_rank() const
  {
    return 1;
  }

  std::size_t value_dimension(std::size_t i) const
  {
    return 3;
  }

  void update()
  {}

};

class myPointer
{
  public:
  //myPointer() : pointer ((new dolfin::HDF5File(MPI_COMM_WORLD,"/u/laureandi/ifumagalli/dcp_test_output/impostoVertCos/sol.hdf5","r"))) {}
  myPointer() : pointer ((new dolfin::HDF5File(MPI_COMM_WORLD,"/u/laureandi/ifumagalli/dcp_test_output/tuttaVel/sol.hdf5","r"))) {}

  std::shared_ptr<dolfin::HDF5File> pointer;
};
class ImposedFromFile : public dolfin::Expression
{
public:
  void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x) const
  {
/*      const std::vector<double>& coords (funSp_->dofmap()->tabulate_all_coordinates(* funSp_->mesh()));
      std::size_t nearest (0);
      for (std::size_t i = 0; i!=coords.size()/2; ++i)
      {
        if ( (x[0]-coords[2*i])*(x[0]-coords[2*i])+(x[1]-coords[2*i+1])*(x[1]-coords[2*i+1])
              <
             (x[0]-coords[2*nearest])*(x[0]-coords[2*nearest])+(x[1]-coords[2*nearest+1])*(x[1]-coords[2*nearest+1]) )
          nearest = i;
      }
      
*/
/*        if (timeStep_.val == 40)
          * fun_->vector() = accumulatedValues_;
        else
          (* fun_->vector()) *= 0;*/
        (* fun_)[0].eval (values,x);
  }
  std::size_t value_rank() const
  {
    return 1;
  }

  std::size_t value_dimension(std::size_t i) const
  {
    return 2;
  }

/*  void setFunSp (dolfin::FunctionSpace& funSp)
  {
      funSp_.reset (&funSp);
  }*/
  void setMesh (const dolfin::Mesh& mesh)
  {
      funSp_.reset (new myNavierstokesTimeCurvLinear::FunctionSpace(mesh));
  }
  void loadFun (std::string funname)
  {
      fun_.reset (new dolfin::Function(*funSp_));
      sourceFile_.pointer->read(*fun_,funname);
  }
  void loadFun ()
  {
      fun_.reset (new dolfin::Function(*funSp_));
      sourceFile_.pointer->read(*fun_,std::to_string(timeStep_.val));
      if (accumulatedValues_.empty())
      {
        accumulatedValues_.init(MPI_COMM_WORLD,fun_->vector()->size());
        accumulatedValues_ = * fun_->vector();
      }
      else
        accumulatedValues_ += * fun_->vector();
  }
  void update ()
  {
    ++(timeStep_.val);
  }

private:
    myPointer sourceFile_;
    std::shared_ptr<dolfin::FunctionSpace> funSp_;
    std::shared_ptr<dolfin::Function> fun_;
    dolfin::Vector accumulatedValues_;
    myInt timeStep_;
};

class FunctionToExpression : public dolfin::Expression
{
public:
  
  explicit FunctionToExpression (std::shared_ptr<dolfin::Function> fun, std::size_t valueDim=1, std::size_t component=0):
      fun_ (fun), leftLength_ (ly), rightLength_ (ly), valueDim_ (valueDim), component_ (component)
      {}
  
  explicit FunctionToExpression (const FunctionToExpression & other) :
      fun_ (other.fun_), leftLength_ (ly), rightLength_ (ly), valueDim_ (other.valueDim_), component_ (other.component_)
      {
  std::cerr << "hello" << fun_->value_size() << std::endl;
      }

  void setLengths (double leftLength, double rightLength)
  {
      leftLength_ = leftLength;
      rightLength_ = rightLength;
  }

  void eval (dolfin::Array<double>& values, const dolfin::Array<double>& x) const
  {
//std::cerr << " ut OK fino a " << __LINE__ << std::endl;
      double xRef ( x[1] / (dolfin::near(x[0],0) ? leftLength_ : rightLength_) );
//std::cerr << " ut OK fino a " << __LINE__ << std::endl;
      dolfin::Array<double> xx (1,&xRef);
      dolfin::Array<double> val (valueDim_);
//std::cerr << " ut OK fino a " << __LINE__ << std::endl;
//std::cerr << fun_->value_size() << ' ' << valueDim_ << ' ' << val.size() << ' ' << xx.size() << std::endl;
      fun_->eval (val,xx);
//std::cerr << " ut OK fino a " << __LINE__ << std::endl;
      for (std::size_t i=0; i!=valueDim_; ++i)
          values[i] = (i==component_ ? val[0] : 0);
  }

  std::size_t value_rank() const
  {
    return valueDim_-1;
  }

  std::size_t value_dimension(std::size_t i) const
  {
    return valueDim_;
  }

  void updateRef (double scale)
  {
      leftLength_ += scale*(*this->fun_)(1);
      rightLength_ = leftLength_;//(*this)(lx,rightLength_);
  }

  double length ()
  {
      return leftLength_;
  }

private:
  //dolfin::Function fun_;
  std::shared_ptr<dolfin::Function> fun_;
  double leftLength_,rightLength_;
  std::size_t valueDim_, component_;
};

class FunctionToExpressionBis : public dolfin::Expression
{
public:
  explicit FunctionToExpressionBis (std::shared_ptr<const dolfin::Function> fun, std::size_t component=0, double length=ly):
      fun_ (fun), component_ (component), length_ (length)
      {}
  
  void eval (dolfin::Array<double>& values, const dolfin::Array<double>& x) const
  {
    double xRef [2] = {0,x[0] * length_};
    dolfin::Array<double> xx (2,xRef);
    (* fun_)[component_].eval (values,xx);
  }

  void updateRef ()
  {
      length_ += (*this)(0,1);
  }

private:
  std::shared_ptr<const dolfin::Function> fun_;
  std::size_t component_;
  double length_;
};

class FirstExtr : public dolfin::SubDomain
{
  public :
    bool inside (const dolfin::Array<double>& x, bool on_boundary) const
    {
      return dolfin::near(x[0],0);
    }
};

class SecondExtr : public dolfin::SubDomain
{
  public :
    bool inside (const dolfin::Array<double>& x, bool on_boundary) const
    {
      return dolfin::near(x[0],1);
    }
};

bool compare (std::tuple<std::size_t,double,double> p1, std::tuple<std::size_t,double,double> p2);
bool compareDofs (std::tuple<std::size_t,double,double,std::size_t> p1, std::tuple<std::size_t,double,double,std::size_t> p2);
bool compareFirst (std::tuple<std::size_t,double,double> p1, std::tuple<std::size_t,double,double> p2);
std::size_t getFirst (const std::tuple<std::size_t,double,double> t);
std::size_t getFirstDofs (const std::tuple<std::size_t,double,double,std::size_t> t);

class NormalAtVertex;

#endif
