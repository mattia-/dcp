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

#include <dolfin.h>

// parameters
#define yxratio 2
//#define lx 0.05
#define nx 20
#define lx 0.05
#define ly lx*yxratio
//#define ustar 0.005
#define ustar 0.02
#define TP_EPS lx/nx+6e-16
/*double yxratio = 2;
double lx (0.05), ly (lx*yxratio);
double ustar (0.005);
double TP_EPS (lx/40+6e-16);*/

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
    bool inside (const dolfin::Array<double>& x, bool on_boundary) const
    {
        return on_boundary
               && ( dolfin::near (x[0],0.0) || dolfin::near(x[0],lx) );
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
              && ( dolfin::near(x[0],0) && dolfin::near(x[1],ly) );
    }
};
class TriplePointRightVertex : public dolfin::SubDomain
{
  public:
    bool inside (const dolfin::Array<double>& x, bool on_boundary) const
    {
        return on_boundary
              && ( dolfin::near(x[0],lx) && dolfin::near(x[1],ly) );
    }
};
class TopBoundary : public dolfin::SubDomain
{
    bool inside (const dolfin::Array<double>& x, bool on_boundary) const
    {
        return on_boundary
               &&
              ( x[1] >= ly - 6.0e-16 );
/*          //what follows takes away the triple points
               &&
             !( dolfin::near(x[0],0) )
               &&
             !( dolfin::near(x[0],lx) );*/
    }
};
class BottomBoundary : public dolfin::SubDomain
{
    bool inside (const dolfin::Array<double>& x, bool on_boundary) const
    {
        return on_boundary
               &&
              (dolfin::near (x[1],0.0));
    }
};
class InflowDirichletData : public dolfin::Expression
{
    void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x) const
    {
        values[0] = 0.0;
        values[1] = ustar;// * 4/(lx*lx)*x[0]*(lx-x[0]) ;
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
        values[2] = 10;
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

#endif
