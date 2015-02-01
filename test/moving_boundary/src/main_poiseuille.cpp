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
#include "poiseuille.h"
#include <differential_problems/differential_problems.h>
#include "MovingTimeDependentProblem.h"
#include "MovingTimeDependentProblem.cpp"
//#include "MovingTimeDependentProblem.h"
//#include "MovingTimeDependentProblem.cpp"

namespace poiseuille
{
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
            poiseuille::Circle circle;
    }; 
*/
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
            return (dolfin::near (x[0], 4) && on_boundary);
        } 
    };
    
    class BottomBoundary : public dolfin::SubDomain
    {
        bool inside (const dolfin::Array<double>& x, bool on_boundary) const
        {
            return (x[1]<0.4 && !(dolfin::near (x[0],0)|| dolfin::near(x[0],4)) && on_boundary);
        }
    };

    class TopBoundary : public dolfin::SubDomain
    {
        bool inside (const dolfin::Array<double>& x, bool on_boundary) const
        {
            return (x[1]>0.6 && !(dolfin::near (x[0],0)|| dolfin::near(x[0],4)) && on_boundary);
/*            return on_boundary
		               && !(
                      dolfin::near(x[0], 0)
           		     ||
             		      dolfin::near(x[0], 1)
          		     ||
              		    dolfin::near(x[1], 0));  */
        }
    };

    class Inflow : public dolfin::Expression
    {
        void eval (dolfin::Array<double>& values, const dolfin::Array<double>& x) const
        {
            values[0] = (1+x[1])*(1-x[1]);
            values[1] = 0;
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

    class InitialSolution : public dolfin::Expression
    {
        void eval (dolfin::Array<double>& values, const dolfin::Array<double>& x) const
        {
            values[0] = values[0] = (1+x[1])*(1-x[1]);
            values[1] = 0;
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
}

int main (int argc, char* argv[])
{

/*dolfin::UnitSquareMesh msh(10,10);
geometry::MapTgamma mapT;
mapT.initTime(0);
mapT.setTime(2);
dolfin::ALE mshMover;
dolfin::Mesh tmpMsh(msh);
dolfin::plot(tmpMsh,"tmp init"); dolfin::interactive();
mshMover.move(tmpMsh,mapT);
dolfin::plot(tmpMsh,"tmp moved"); dolfin::interactive();
dolfin::plot(msh,"mesh init"); dolfin::interactive();
std::shared_ptr<dolfin::MeshDisplacement> disp (mshMover.move(msh,tmpMsh));
dolfin::plot(msh,"mesh moved"); dolfin::interactive();
dolfin::plot(tmpMsh,"tmp moved moved"); dolfin::interactive();

if (true)
  return 0;*/

    // create mesh and finite element space
    std::cout << "Create mesh and finite element space..." << std::endl;
/*    mshr::Rectangle rectangle (dolfin::Point (0.0, 0.0), dolfin::Point (10.0, 7.0));
    mshr::Circle circle (dolfin::Point (3.5, 3.5), 0.5);
    dolfin::Mesh mesh;
    mshr::generate (mesh, *(rectangle - circle), 50);
*/
    dolfin::RectangleMesh mesh (0,0, 4,1, 40,10);

/*    dolfin::Mesh mesh2(mesh);
    geometry::MapTgamma mapT;
    mapT.initTime(0.0);
    dolfin::ALE::move(mesh,mapT);
    dolfin::plot(mesh);
    mapT.setTime(0.7);
    dolfin::ALE::move(mesh,mapT);
    dolfin::plot(mesh);
    mapT.setTime(0.0); dolfin::ALE::move(mesh,mapT);
    geometry::MeshManager<dolfin::ALE,dolfin::FunctionSpace> meshMan(std::shared_ptr<dolfin::ALE>(new dolfin::ALE()),std::shared_ptr<dolfin::FunctionSpace>(new poiseuille::FunctionSpace(mesh2)));
    mapT.initTime(0.0);
    mapT.setTime(0.3);
    meshMan.moveMesh(mapT);
    mapT.printTimes();
    dolfin::plot(* meshMan.functionSpace().mesh());
    mapT.setTime(0.7);
    meshMan.moveMesh(mapT);
    mapT.printTimes();
    dolfin::plot(* meshMan.functionSpace().mesh());
    mapT.setTime(1.4);
    meshMan.moveMesh(mapT);
    mapT.printTimes();
    dolfin::plot(* meshMan.functionSpace().mesh());
    dolfin::interactive();
    
if (false)
  return 0;*/

    poiseuille::FunctionSpace V (mesh);
    
    // define problem
    double t0 = 0.0;
    double dt = 0.1;
    double T = 3.0;
    std::cout << "Define the problem..." << std::endl;
    Ivan::MovingTimeDependentProblem navierStokesProblem (dolfin::reference_to_no_delete_pointer (mesh),
                                                   dolfin::reference_to_no_delete_pointer (V),
                                                   t0,
                                                   dt, 
                                                   T, 
                                                   std::vector<std::string> ({"residual_form", "jacobian_form"}),
                                                   std::vector<std::string> ({"residual_form"})
                                                  );
         
    dcp::NonlinearProblem <poiseuille::ResidualForm, poiseuille::JacobianForm> 
        timeSteppingProblem (dolfin::reference_to_no_delete_pointer (mesh), 
                             dolfin::reference_to_no_delete_pointer (V),
                             "trial");
    navierStokesProblem.setTimeSteppingProblem (dolfin::reference_to_no_delete_pointer (timeSteppingProblem));
    
    // define constant
    std::cout << "Define the problem's coefficients..." << std::endl;
    dolfin::Constant nu (1e-1);
//    dolfin::Constant inflowDirichletBC (1.0, 0.0);
    poiseuille::Inflow inflowDirichletBC;
    dolfin::Constant symmetryDirichletBC (0.0);
    dolfin::Constant noSlipDirichletBC (0.0, 0.0);

    // define dirichlet boundary conditions
    std::cout << "Define the problem's Dirichlet boundary conditions..." << std::endl;
/*    poiseuille::InflowBoundary inflowBoundary;
    poiseuille::GammaSD gammaSD;
    poiseuille::NoSlipBoundary noSlipBoundary;
    
    navierStokesProblem.addDirichletBC (dolfin::DirichletBC (*V[0], inflowDirichletBC, inflowBoundary));
    navierStokesProblem.addDirichletBC (dolfin::DirichletBC (*V[0], noSlipDirichletBC, noSlipBoundary, "topological", false));
    navierStokesProblem.addDirichletBC (dolfin::DirichletBC (*(*V[0])[1], symmetryDirichletBC, gammaSD));
*/
    poiseuille::LeftBoundary inflowBoundary;
    poiseuille::RightBoundary gammaSD;
    poiseuille::TopBoundary noSlipBoundaryTop;
    poiseuille::BottomBoundary noSlipBoundaryBottom;
    navierStokesProblem.addDirichletBC (dolfin::DirichletBC (*V[0], inflowDirichletBC, inflowBoundary));
    navierStokesProblem.addDirichletBC (dolfin::DirichletBC (*V[0], noSlipDirichletBC, noSlipBoundaryTop, "topological", false));
    navierStokesProblem.addDirichletBC (dolfin::DirichletBC (*(*V[0])[1], symmetryDirichletBC, noSlipBoundaryBottom, "topological", false));
//    navierStokesProblem.addDirichletBC (dolfin::DirichletBC (*(*V[0])[1], symmetryDirichletBC, gammaSD));

    // problem settings
    std::cout << "Set the problem's coefficients..." << std::endl;
    navierStokesProblem.setCoefficient ("residual_form", dolfin::reference_to_no_delete_pointer (nu), "nu");
    navierStokesProblem.setCoefficient ("jacobian_form", dolfin::reference_to_no_delete_pointer (nu), "nu");
    
    navierStokesProblem.parameters ["time_stepping_solution_component"] = 0;
    
    //navierStokesProblem.setInitialSolution (dolfin::Constant (1,0,0));
    navierStokesProblem.setInitialSolution (poiseuille::InitialSolution());

    // solve problem
    std::cout << "Solve the problem..." << std::endl;
    navierStokesProblem.solve ();

    // plots
    dolfin::plot (* navierStokesProblem.solution().function_space()->mesh(), "Final mesh");
    
    dolfin::interactive ();
    
    return 0;
}
