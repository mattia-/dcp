/* 
 *  Copyright (C) 2014, Mattia Tamellini, mattia.tamellini@gmail.com
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

#include <iostream>
#include <string>
#include <dolfin.h>
//#include <mshr.h>
#include <differential_problems/differential_problems.h>
#include "navierstokes.h"
#include "geometry.h"
#include "IvanNonlinearProblem.h"

//bool mydebug = false;

namespace navierstokes
{
    class LeftBoundary : public dolfin::SubDomain
    {
        bool inside (const dolfin::Array<double>& x, bool on_boundary) const
        {
            return (dolfin::near (x[0], 0) && on_boundary);
        } 
    };
    
    class RightBoundary : public dolfin::SubDomain
    {
        bool inside (const dolfin::Array<double>& x, bool on_boundary) const
        {
            return (dolfin::near (x[0], 1) && on_boundary);
        } 
    };
    
    class BottomBoundary : public dolfin::SubDomain
    {
        bool inside (const dolfin::Array<double>& x, bool on_boundary) const
        {
            return dolfin::near (x[1], 0) && on_boundary;
        }
    };

    class TopBoundary : public dolfin::SubDomain
    {
        bool inside (const dolfin::Array<double>& x, bool on_boundary) const
        {
            return on_boundary
		               && !(
                      dolfin::near(x[0], 0)
           		     ||
             		      dolfin::near(x[0], 1)
          		     ||
              		    dolfin::near(x[1], 0));
        }
    };

    class SourceTerm : public dolfin::Expression
    {
        void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x) const
        {
          values[0] = 0;//1000000.0*(x[1]>1);
        }
/*        uint value_rank() const
        {
          return 1;
        }
        uint value_dimension(uint i) const
        {
          return 2;
        }
*/  };

/*    class Inflow : public dolfin::Expression
    {
        void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x) const
        {
          values[0] = 0;
          values[1] = 1-x[1];
        }
        std::size_t value_rank() const
        {
          return 1;
        }
        std::size_t value_dimension(uint i) const
        {
          return 2;
        }
        
    };
*/
}
/*    class InflowBoundary : public dolfin::SubDomain
    {
        bool inside (const dolfin::Array<double>& x, bool on_boundary) const
        {
            return x[0] <= (0 + DOLFIN_EPS) && on_boundary;
        }
    };

    class GammaSD : public dolfin::SubDomain
    {
        bool inside (const dolfin::Array<double>& x, bool on_boundary) const
        {
            return (x[1] <= (0 + DOLFIN_EPS) || x[1] >= (7 - DOLFIN_EPS)) && on_boundary;
        }
    };
    
    class Circle : public dolfin::SubDomain
    {
        public:
            Circle () : 
                center (2),
                radius (0.5)
            {
                center[0] = 3.5;
                center[1] = 3.5;
            }
            
            bool inside (const dolfin::Array<double>& x, bool on_boundary) const
            {
                double dx = x[0] - center[0];
                double dy = x[1] - center[1];
                double r = sqrt (dx * dx + dy * dy);
                
                return r <= (radius + 1e-3) && on_boundary;
            }

        private:
            dolfin::Array<double> center;
            double radius;
    }; 
    
    class NoSlipBoundary : public dolfin::SubDomain
    { 
        public:
            bool inside (const dolfin::Array<double>& x, bool on_boundary) const
            {
                bool onWall1 = (x[0] >= 1.75 - DOLFIN_EPS && x[0] <= 1.75 + DOLFIN_EPS)
                            && (x[1] >= 3 - DOLFIN_EPS && x[1] <= 4 + DOLFIN_EPS)
                            && on_boundary;
                bool onWall2 = (x[0] >= (1.75 - DOLFIN_EPS) && x[0] <= (2 + DOLFIN_EPS))
                            && (x[1] <= (4 + DOLFIN_EPS) && x[1] >= (4 - DOLFIN_EPS))
                            && on_boundary;
                bool onWall3 = (x[0] >= 2 - DOLFIN_EPS && x[0] <= 2 + DOLFIN_EPS)
                            && (x[1] >= 3 - DOLFIN_EPS && x[1] <= 4 + DOLFIN_EPS)
                            && on_boundary;
                bool onWall4 = (x[0] >= (1.75 - DOLFIN_EPS) && x[0] <= (2 + DOLFIN_EPS))
                            && (x[1] <= (3 + DOLFIN_EPS) && x[1] >= (3 - DOLFIN_EPS))
                            && on_boundary;
                
                bool onCircle = circle.inside (x, on_boundary);

                return onWall1 || onWall2 || onWall3 || onWall4 || onCircle;
            }

        private:
            navierstokes::Circle circle;
    }; 
}*/

int main (int argc, char* argv[])
{

if(mydebug==true){
dolfin::set_log_level(dolfin::DBG);
}
    // create mesh and finite element space
    std::cout << "Create mesh and finite element space..." << std::endl;
/*    mshr::Rectangle rectangle (dolfin::Point (0.0, 0.0), dolfin::Point (10.0, 7.0));
    mshr::Circle circle (dolfin::Point (3.5, 3.5), 0.5);
    dolfin::Mesh mesh;
    mshr::generate (mesh, *(rectangle - circle), 50);
*/
    dolfin::UnitSquareMesh mesh (10,10);

    // define map moving the domain
    geometry::MapTgamma T;
    T.init_count();

    dolfin::File solutionFile("solution.pvd");
    
    navierstokes::FunctionSpace Vold (mesh);
    
    // define problem
    std::cout << "Define the problem..." << std::endl;
    Ivan::NonlinearProblem <navierstokes::ResidualForm, navierstokes::JacobianForm> 
//    dcp::NonlinearProblem <navierstokes::ResidualForm, navierstokes::JacobianForm> 
        navierStokesProblem (dolfin::reference_to_no_delete_pointer (mesh), 
                             dolfin::reference_to_no_delete_pointer (Vold),
                             "trial");

    // the mesh mover
    geometry::MeshManager<dolfin::ALE,dolfin::FunctionSpace> meshManager(std::shared_ptr<dolfin::ALE>(new dolfin::ALE()),navierStokesProblem.functionSpace());

for(std::size_t i=1; i<=2; i++)
{
    
    const dolfin::FunctionSpace& V (meshManager.functionSpace());

    // define constant
    std::cout << "Define the problem's coefficients..." << std::endl;
    dolfin::Constant nu (1e-1);
    dolfin::Constant inflowDirichletBC (1.0, 0.0);
//    navierstokes::Inflow inflowDirichletBC;
    dolfin::Constant symmetryDirichletBC (0.0);
    dolfin::Constant noSlipDirichletBC (0.0, 0.0);

    // define dirichlet boundary conditions
    std::cout << "Define the problem's Dirichlet boundary conditions..." << std::endl;
    navierstokes::LeftBoundary inflowBoundary;
    navierstokes::RightBoundary gammaSD;
    navierstokes::TopBoundary noSlipBoundaryTop;
    navierstokes::BottomBoundary noSlipBoundaryBottom;
    
    navierStokesProblem.addDirichletBC (dolfin::DirichletBC (*V[0], inflowDirichletBC, inflowBoundary));
    navierStokesProblem.addDirichletBC (dolfin::DirichletBC (*V[0], noSlipDirichletBC, noSlipBoundaryTop, "topological", false));
    navierStokesProblem.addDirichletBC (dolfin::DirichletBC (*V[0], noSlipDirichletBC, noSlipBoundaryBottom, "topological", false));
    navierStokesProblem.addDirichletBC (dolfin::DirichletBC (*(*V[0])[1], symmetryDirichletBC, gammaSD));

    // problem settings
    std::cout << "Set the problem's coefficients..." << std::endl;
    navierStokesProblem.setCoefficient ("residual_form", dolfin::reference_to_no_delete_pointer (nu), "nu");
    navierStokesProblem.setCoefficient ("jacobian_form", dolfin::reference_to_no_delete_pointer (nu), "nu");

if (mydebug==true){
std::cout << dolfin::get_log_level() << '\t';
//dolfin::set_log_level(dolfin::DBG);
std::cout << dolfin::get_log_level() << std::endl;
}

    // solve problem
    std::cout << "Solve the problem..." << std::endl;
    navierStokesProblem.solve ();
if(mydebug==true){
std::cerr << " line : " << __LINE__ << std::endl;
}

    // plots
if(mydebug==true){
    dolfin::plot (mesh, "Mesh");
std::cerr << " line : " << __LINE__ << std::endl;
}
    dolfin::plot (navierStokesProblem.solution ()[0], "Velocity");
if(mydebug==true){
std::cerr << " line : " << __LINE__ << std::endl;
}
    dolfin::plot (navierStokesProblem.solution ()[1], "Pressure");
if(mydebug==true){
std::cerr << " line : " << __LINE__ << std::endl;
}
    dolfin::interactive ();
    
    // move mesh
    T.update_count();
    meshManager.moveMesh(T);
//    meshMover.move(*(navierStokesProblem.mesh()),T);
if(mydebug==true){
    dolfin::plot(navierStokesProblem.mesh(),"dal problem - reference");
    dolfin::interactive();
    dolfin::plot(navierStokesProblem.functionSpace()->mesh(),"dal problem functionSpace - reference");
    dolfin::interactive();
    dolfin::plot(navierStokesProblem.solution().function_space()->mesh(),"da problem solution functionSpace - reference");
    dolfin::interactive();
/*    dolfin::plot(*(navierStokesProblem.mesh()), "dal problem - pointer");
    dolfin::interactive();
*/// col pointer non c'è più
    dolfin::plot(*(V.mesh()),"dal function space");
    dolfin::interactive();
 }

//    V = navierStokes::FunctionSpace(*navierStokesProblem.mesh());
/*    dolfin::plot(*(V.mesh()),"dal function space - dopo =");
    dolfin::interactive();
//    navierStokesProblem.setFunctionSpace(dolfin::reference_to_no_delete_pointer(V));
    dolfin::plot(navierStokesProblem.mesh(),"dal problem - reference - dopo =");
    dolfin::interactive();
*/
/*    dolfin::plot(*(navierStokesProblem.mesh()), "dal problem - pointer - dopo =");
    dolfin::interactive();
*/// col pointer non c'è più

}

    return 0;
}
