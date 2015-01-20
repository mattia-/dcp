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
#include <DifferentialProblem/DifferentialProblem.hpp>
#include "poisson.h"
#include "geometry.h"
#include "MovingFunctionSpace_tmp.h"
#include "IvanLinearDifferentialProblem.hpp"

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
}

int main (int argc, char* argv[])
{
    // create mesh and finite element space 
    std::cout << "Create mesh and finite element space..." << std::endl;
    dolfin::UnitSquareMesh meshold (10, 10);
    poisson::FunctionSpace Vold (meshold);

    // define map moving the domain
    geometry::MapTgamma T;
    T.init_count();
    T.update_count();

    // IL MIO FUNCTIONSPACE
//    geometry::FunctionSpace V(oldV);
    dolfin::ALE meshMover;

    //std::shared_ptr<const dolfin::FiniteElement> element(V.element());
    geometry::FunctionSpace<poisson::FunctionSpace, dolfin::ALE> V (Vold,meshMover);
//    geometry::FunctionSpace<poisson::FunctionSpace, dolfin::ALE> V2 (Vold,meshMover);
/*    std::shared_ptr<geometry::FunctionSpace<poisson::FunctionSpace, dolfin::ALE> > pMfs (dolfin::reference_to_no_delete_pointer(V));
    std::shared_ptr<dolfin::FunctionSpace> pDfs (dolfin::reference_to_no_delete_pointer(V));
    std::shared_ptr<geometry::FunctionSpace<poisson::FunctionSpace, dolfin::ALE> > pMfs_from_pDfs
        (static_cast<std::shared_ptr<geometry::FunctionSpace<poisson::FunctionSpace, dolfin::ALE> > > (pDfs);
*/  geometry::FunctionSpace<poisson::FunctionSpace, dolfin::ALE>*  pMfs (&V);
    dolfin::FunctionSpace* pDfs (&V);
    geometry::FunctionSpace<poisson::FunctionSpace, dolfin::ALE>* pMfs_from_pDfs
        (static_cast<geometry::FunctionSpace<poisson::FunctionSpace, dolfin::ALE>* > (pDfs));
    dolfin::plot(V.mesh(),"function space 1");
    dolfin::interactive();
    dolfin::plot(pMfs->mesh(),"ptr mio 1");
    dolfin::interactive();
    dolfin::plot(pDfs->mesh(),"ptr dolfin 1");
    dolfin::interactive();
    dolfin::plot(pMfs_from_pDfs->mesh(),"ptr mio da dolfin 1");
    dolfin::interactive();
    V.moveMesh(T);
    dolfin::plot(V.mesh(),"function space 2");
    dolfin::interactive();
    dolfin::plot(pMfs->mesh(),"ptr mio 2");
    dolfin::interactive();
    dolfin::plot(pDfs->mesh(),"ptr dolfin 2");
    dolfin::interactive();
    dolfin::plot(pMfs_from_pDfs->mesh(),"ptr mio da dolfin 2");
    dolfin::interactive();

if(false)
return 0;

    dolfin::File solutionFile("solution.pvd");

for(std::size_t i=1; i<=4; i++)
{
    // define problem
    std::cout << "Define the problem..." << std::endl;
    dcp::IvanLinearDifferentialProblem <poisson::BilinearForm, poisson::LinearForm> 
        poissonProblem (V.mesh(), 
//        poissonProblem (dolfin::reference_to_no_delete_pointer (mesh), 
                        dolfin::reference_to_no_delete_pointer (V));

    // define coefficients
    std::cout << "Define the problem's coefficients..." << std::endl;
    dolfin::Constant k (1.0);
    dolfin::Constant f (1.0);
    dolfin::Constant g (1.0);
    dolfin::Constant hbot (0.0);
    dolfin::Constant htop (10.0);

    // define dirichlet boundary conditions 
    std::cout << "Define the problem's Dirichlet boundary conditions..." << std::endl;
    poisson::BottomBoundary bottomEdge;
    poissonProblem.addDirichletBC (dolfin::DirichletBC (V, hbot, bottomEdge));
    poisson::TopBoundary topEdge;
    poissonProblem.addDirichletBC (dolfin::DirichletBC (V, htop, topEdge));
    
    // define neumann boundary conditions 
    std::cout << "Define the problem's Neumann boundary conditions..." << std::endl;
    poisson::LateralBoundary neumannBoundary;
    dolfin::FacetFunction<std::size_t> meshFacets (*(V.mesh()));
    meshFacets.set_all (0);
    neumannBoundary.mark (meshFacets, 1);
    
    std::cout << dolfin::get_log_level() << '\t';
//    dolfin::set_log_level(dolfin::DBG);
    std::cout << dolfin::get_log_level() << std::endl;

    // set problem coefficients
    std::cout << "Set the problem's coefficients..." << std::endl;
    poissonProblem.setCoefficient ("bilinear_form", dolfin::reference_to_no_delete_pointer (k), "k");
    poissonProblem.setCoefficient ("linear_form", dolfin::reference_to_no_delete_pointer (k), "k");
    poissonProblem.setCoefficient ("linear_form", dolfin::reference_to_no_delete_pointer (f), "f");
    poissonProblem.setCoefficient ("linear_form", dolfin::reference_to_no_delete_pointer (g), "g");

    // solve problem
    std::cout << "Solve the problem..." << std::endl;
    poissonProblem.solve (true);

    // plots
    dolfin::Function const u = poissonProblem.solution();
    solutionFile << u;
    dolfin::plot (u, "Solution");
/*    dolfin::VTKPlotter plotter(std::shared_ptr<const dolfin::Function>(&u));
//    dolfin::plot (mesh);
    plotter.elevate (0.0);
    plotter.plot ();
*/
/*    dolfin::plot (mesh, "Mesh");
    auto vtk_plotter = dolfin::plot (poissonProblem.solution(), "Solution");
    vtk_plotter -> elevate(0.0);
    vtk_plotter -> plot ();//poissonProblem.solution());
    
    dolfin::interactive ();*/

    // verifiche del caso    
    std::cout << "function space : " << V.str(false) << "\nmesh : " << V.mesh()->str(false) << std::endl;
//    std::shared_ptr<const dolfin::Mesh> meshnew(Vnew.mesh());
//    auto meshMovernew(Vnew.meshMover());

    // move mesh
    T.update_count();
//    meshMovernew->move(mesh,T);
//    meshMover.move(mesh,T);
//    meshMovernew->move(*meshnew,T);  // non va: meshnew e' const...*
    V.moveMesh(T);
    dolfin::plot(*(V.mesh()),"dal function space");
    dolfin::interactive();
    dolfin::plot(poissonProblem.mesh(),"dal problem");
    dolfin::plot(u.function_space()->mesh(),"dalla solution");
    dolfin::interactive();
    const dolfin::Mesh& meshCastTmp (*(u.function_space()->mesh()));
    dolfin::Mesh& meshCast (const_cast<dolfin::Mesh&> (meshCastTmp));
    dolfin::plot(meshCast,"dal cast");
    dolfin::interactive();
    meshMover.move(meshCast,T);
    dolfin::plot(meshCast,"dal cast - moved");
    dolfin::interactive();

//    V = poisson::FunctionSpace(mesh);
//    poisson::FunctionSpace V2(meshnew);
//    dolfin::plot (mesh);
//    std::cin.get();
}

/*    std::cout << "function space : " << V2.str(false) << "\nmesh : " << V2.mesh()->str(false) << std::endl;
    dolfin::plot(*(V2.mesh()),"original V2 mesh");
    V2 = V;
    std::cout << "function space : " << V2.str(false) << "\nmesh : " << V2.mesh()->str(false) << std::endl;
    dolfin::plot(*(V2.mesh()),"assigned V2 mesh");
    dolfin::interactive();
*/

    return 0;
}
