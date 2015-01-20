/* 
 *  Copyright (C) 2015, Ivan Fumagalli, ivanfumagalli.if@gmail.com
 * 
 *  Modified from linear_differential_problems/main.cpp
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
#include <differential_problems/differential_problems.h>
#include "poisson.h"
#include "geometry.h"

bool mydebug = false;

namespace poisson
{
    class LateralBoundary : public dolfin::SubDomain
    {
        bool inside (const dolfin::Array<double>& x, bool on_boundary) const
        {
            return (dolfin::near (x[0], 0) && on_boundary)
                   ||
                   (dolfin::near (x[0], 1) && on_boundary);
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
}

int main (int argc, char* argv[])
{
    // create mesh and finite element space 
    std::cout << "Create mesh and finite element space..." << std::endl;
    dolfin::UnitSquareMesh mesh (100, 100);
    poisson::FunctionSpace Vp (mesh);

    // define map moving the domain
    geometry::MapTgamma T;
    T.init_count();
//    T.update_count();

    // the mesh mover
//    dolfin::ALE meshMover;
    dolfin::File solutionFile("solution.pvd");

    // define problem
    std::cout << "Define the problem..." << std::endl;
    dcp::LinearProblem <poisson::BilinearForm, poisson::LinearForm> 
        poissonProblem (dolfin::reference_to_no_delete_pointer (mesh), 
                        dolfin::reference_to_no_delete_pointer (Vp));

if (mydebug==true){
std::cerr << "linea " <<  __LINE__ << std::endl;
}
/*    // define mesh manager
    std::shared_ptr<poisson::FunctionSpace> ptr (dolfin::reference_to_no_delete_pointer(static_cast<poisson::FunctionSpace&>(*poissonProblem.functionSpace())));
if (mydebug==true){
std::cerr << "linea " <<  __LINE__ << std::endl;
}
dolfin::plot(*ptr->mesh());
if (mydebug==true){
std::cerr << "linea " <<  __LINE__ << std::endl;
}
    geometry::MeshManager<dolfin::ALE,poisson::FunctionSpace> meshManager(std::shared_ptr<dolfin::ALE>(&meshMover),ptr);
*/    geometry::MeshManager<dolfin::ALE,dolfin::FunctionSpace> meshManager(std::shared_ptr<dolfin::ALE>(new dolfin::ALE()),poissonProblem.functionSpace());
  // ??? la riga seguente dà segmentation fault al costruttore di meshMover_ del meshManager: perché???
    //geometry::MeshManager<dolfin::ALE,dolfin::FunctionSpace> meshManager(std::shared_ptr<dolfin::ALE>(&meshMover),poissonProblem.functionSpace());
if (mydebug==true){
std::cerr << "linea " <<  __LINE__ << std::endl;
}

    const dolfin::FunctionSpace& V (meshManager.functionSpace());

    // define coefficients
    std::cout << "Define the problem's coefficients..." << std::endl;
    dolfin::Constant k (1.0);
//    dolfin::Constant f (1.0);
    poisson::SourceTerm f;
    dolfin::Constant g (1.0);
    dolfin::Constant hbot (0.0);
    dolfin::Constant htop (10.0);

    // define dirichlet boundary conditions 
    std::cout << "Define the problem's Dirichlet boundary conditions..." << std::endl;
    poisson::BottomBoundary bottomEdge;
if (mydebug==true){
std::cerr << "linea " <<  __LINE__ << std::endl;
}
    poissonProblem.addDirichletBC (dolfin::DirichletBC (V, hbot, bottomEdge));
if (mydebug==true){
std::cerr << "linea " <<  __LINE__ << std::endl;
}
    poisson::TopBoundary topEdge;
    poissonProblem.addDirichletBC (dolfin::DirichletBC (V, htop, topEdge));
if (mydebug==true){
std::cerr << "linea " <<  __LINE__ << std::endl;
}
    
    // define neumann boundary conditions 
    std::cout << "Define the problem's Neumann boundary conditions..." << std::endl;
    poisson::LateralBoundary neumannBoundary;
if (mydebug==true){
std::cerr << "linea " <<  __LINE__ << std::endl;
}
    dolfin::FacetFunction<std::size_t> meshFacets (*(V.mesh()));
if (mydebug==true){
std::cerr << "linea " <<  __LINE__ << std::endl;
}
    meshFacets.set_all (0);
if (mydebug==true){
std::cerr << "linea " <<  __LINE__ << std::endl;
}
    neumannBoundary.mark (meshFacets, 1);
if (mydebug==true){
std::cerr << "linea " <<  __LINE__ << std::endl;
}
    
if (mydebug==true){
std::cout << dolfin::get_log_level() << '\t';
//dolfin::set_log_level(dolfin::DBG);
std::cout << dolfin::get_log_level() << std::endl;
}

    // set problem coefficients
    std::cout << "Set the problem's coefficients..." << std::endl;
    poissonProblem.setCoefficient ("bilinear_form", dolfin::reference_to_no_delete_pointer (k), "k");
    poissonProblem.setCoefficient ("linear_form", dolfin::reference_to_no_delete_pointer (k), "k");
    poissonProblem.setCoefficient ("linear_form", dolfin::reference_to_no_delete_pointer (f), "f");
    poissonProblem.setCoefficient ("linear_form", dolfin::reference_to_no_delete_pointer (g), "g");
if (mydebug==true){
std::cerr << "linea " <<  __LINE__ << std::endl;
}

for(std::size_t i=1; i<=5; i++)
{
//    const dolfin::FunctionSpace& V (meshManager.functionSpace());

    // solve problem
    std::cout << "Solve the problem..." << std::endl;
    poissonProblem.solve ();
if (mydebug==true){
dolfin::info(poissonProblem.linearOperator(),true);
}
    solutionFile << poissonProblem.solution();
    
    // plots
    dolfin::plot (mesh);
    dolfin::plot (poissonProblem.solution ());
    
    dolfin::interactive ();
    
if (mydebug==true){
std::cerr << "linea " <<  __LINE__ << std::endl;
}
    // move mesh
    T.update_count();
    meshManager.moveMesh(T);
//    meshMover.move(*(poissonProblem.mesh()),T);
    dolfin::plot(poissonProblem.mesh(),"dal problem - reference");
    dolfin::interactive();
/*    dolfin::plot(*(poissonProblem.mesh()), "dal problem - pointer");
    dolfin::interactive();
*/// col pointer non c'è più
    dolfin::plot(*(V.mesh()),"dal function space");
    dolfin::interactive();

//    V = poisson::FunctionSpace(*poissonProblem.mesh());
    dolfin::plot(*(V.mesh()),"dal function space - dopo =");
    dolfin::interactive();
//    poissonProblem.setFunctionSpace(dolfin::reference_to_no_delete_pointer(V));
    dolfin::plot(poissonProblem.mesh(),"dal problem - reference - dopo =");
    dolfin::interactive();
/*    dolfin::plot(*(poissonProblem.mesh()), "dal problem - pointer - dopo =");
    dolfin::interactive();
*/// col pointer non c'è più

}

    return 0;
}
