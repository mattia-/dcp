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
#include "myNavierstokesTimeCurv.h"
#include <dcp/differential_problems/differential_problems.h>
#include "geometry.h"
#include "MovingTimeDependentProblem.h"
//#include "MovingTimeDependentProblem.cpp"

bool geometry::timeNOTcount(true);

namespace myNavierstokesTimeCurv
{
    // parameters
    double yxratio = 2;
    double lx (0.05), ly (lx*yxratio);
    std::size_t nx (10), ny (nx*yxratio);
    double ustar (0.02);//(0.000001);
    double kinematic_viscosity (1e-6);

    class LeftBoundary : public dolfin::SubDomain
    {
        bool inside (const dolfin::Array<double>& x, bool on_boundary) const
        {
            return dolfin::near (x[0], 0) && on_boundary;
//            return dolfin::near (x[0], 0) && (x[1]<ly-6e-16) && on_boundary;
        } 
    };
    
    class RightBoundary : public dolfin::SubDomain
    {
        bool inside (const dolfin::Array<double>& x, bool on_boundary) const
        {
            return (dolfin::near (x[0], lx) && on_boundary);
//            return dolfin::near (x[0], lx) && (x[1]<ly-6e-16) && on_boundary;
        } 
    };
    
 /*   class BottomBoundary : public dolfin::SubDomain
    {
        bool inside (const dolfin::Array<double>& x, bool on_boundary) const
        {
            return dolfin::near (x[1], 0) && on_boundary;
        }
    };
 */   class BottomBoundaryEvaluator
    { 
        public:
            bool operator() (const dolfin::Array<double>& x, bool on_boundary)
            {
                return dolfin::near (x[1], 0) && on_boundary;
            }
    }; 

    class TopBoundary : public dolfin::SubDomain
    {
        bool inside (const dolfin::Array<double>& x, bool on_boundary) const
        {
            return on_boundary
                   &&
                  (x[1] >= 0.10 - 6.0e-16);
        }
    };

//sistemare ++
	class TriplePoints : public dolfin::SubDomain
    {
        bool inside (const dolfin::Array<double>& x, bool on_boundary) const
        {
            return  on_boundary
                  && ( ( dolfin::near(x[0], 0) && dolfin::near(x[1],ly) )
           		       ||
                       ( dolfin::near(x[0], lx) && dolfin::near(x[1],ly)  ) );
        }
    };
	
/*    class InitialSolution : public dolfin::Expression
    {
        void eval (dolfin::Array<double>& values, const dolfin::Array<double>& x) const
        {
            values[0] = 
        }
    };*/

    class InflowDirichletBC : public dcp::TimeDependentExpression
    {
        public:
            InflowDirichletBC () : dcp::TimeDependentExpression (2) { }
        
            void eval (dolfin::Array<double>& values, const dolfin::Array<double>& x) const override
            {
                values[0] = 0;
                values[1] = sin(2*3.14*t); //x[0]*(1-x[0]);
            }
    }; 

    class InflowDirichletBCEvaluator
    {
        public:
            void operator() (dolfin::Array<double>& values, const dolfin::Array<double>& x, const double& t)
            {
                values[0] = 0;
//                values[1] = sin(2*3.14*t) * 6*x[0]*(1-x[0]); //firstDrop
//                values[1] = 0.5 * sin(2*3.14*t) * 4*x[0]*(1-x[0]); //freq1
                values[1] = ustar * sin(2*3.14*t) * 4/(lx*lx)*x[0]*(lx-x[0]);
            }
    };
    
    class InitialDisplacement : public dolfin::Expression
    {
        void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x) const
        {
            values[0] = 0.0;
            values[1] = 0.0;//05 * x[1] * ( 0.1 * 4/(lx*lx)*x[0]*(lx-x[0]) );
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
    //dolfin::UnitSquareMesh mesh (20,20);
    dolfin::RectangleMesh mesh ( 0.0, 0.0,
                                 myNavierstokesTimeCurv::lx,  myNavierstokesTimeCurv::ly,
                                 myNavierstokesTimeCurv::nx,
                                 myNavierstokesTimeCurv::ny);

    dolfin::ALE meshInitializer;
    myNavierstokesTimeCurv::InitialDisplacement initialDisplacement;
    meshInitializer.move(mesh, initialDisplacement);
    
    myNavierstokesTimeCurv::FunctionSpace V (mesh);

/*dolfin::plot(mesh);
dolfin::plot(* (* V[0])[0]->mesh());
dolfin::interactive();
return 0;
*/
    
    dcp::NonlinearProblem <myNavierstokesTimeCurv::ResidualForm, myNavierstokesTimeCurv::JacobianForm> 
         timeSteppingProblem (dolfin::reference_to_no_delete_pointer (V), 
                             "trial");
    
    // define problem
    double t0 = 0.0;
    double dt = 0.05;
//    double T = 1; //firstDrop
    double T = 2;
    std::cout << "Define the problem..." << std::endl;
 /*   Ivan::MovingTimeDependentProblem navierStokesProblem (dolfin::reference_to_no_delete_pointer (mesh),
                                                   dolfin::reference_to_no_delete_pointer (V),
 */   Ivan::MovingTimeDependentProblem navierStokesProblem (dolfin::reference_to_no_delete_pointer (timeSteppingProblem),
                                                   t0,
                                                   dt, 
                                                   T, 
                                                   {"residual_form", "jacobian_form"},
                                                   {"residual_form"}
                                                  );
         
 /*   dcp::NonlinearProblem <myNavierstokesTimeCurv::ResidualForm, myNavierstokesTimeCurv::JacobianForm> 
  //       timeSteppingProblem (dolfin::reference_to_no_delete_pointer (mesh), 
//treno+        timeSteppingProblem (dolfin::reference_to_no_delete_pointer(const_cast<dolfin::Mesh&> (* navierStokesProblem.functionSpace()->mesh())), 
                             // tutto sto casino serve perché in dolfin::FunctionSpace::mesh() restituisce un std::shared_ptr<const dolfin::Mesh> e il const della Mesh non lo voglio
         timeSteppingProblem (dolfin::reference_to_no_delete_pointer (V), 
  //                          dolfin::reference_to_no_delete_pointer (V),
//treno+                             dolfin::reference_to_no_delete_pointer(* navierStokesProblem.functionSpace()),
//                             navierStokesProblem.functionSpace(),
                             "trial");
    navierStokesProblem.setTimeSteppingProblem (dolfin::reference_to_no_delete_pointer (timeSteppingProblem));
 */
    
    // define constant
    std::cout << "Define the problem's coefficients..." << std::endl;
    dolfin::Constant nu (myNavierstokesTimeCurv::kinematic_viscosity / ( myNavierstokesTimeCurv::ustar * myNavierstokesTimeCurv::lx ));
    dolfin::Constant gamma (7.3e-5);
//    dolfin::Constant p_out (1000.0);
//    dolfin::Constant p_left (0.0);
 //   dolfin::Constant inflowDirichletBC (0.0, 1.0);
    dolfin::Constant symmetryDirichletBC (0.0);
    dolfin::Constant noSlipDirichletBC (0.0, 0.0);

    // define dirichlet boundary conditions
    std::cout << "Define the problem's Dirichlet boundary conditions..." << std::endl;
/*    myNavierstokesTimeCurv::InflowBoundary inflowBoundary;
    myNavierstokesTimeCurv::GammaSD gammaSD;
    myNavierstokesTimeCurv::NoSlipBoundary noSlipBoundary;
    
    navierStokesProblem.addDirichletBC (dolfin::DirichletBC (*V[0], inflowDirichletBC, inflowBoundary));
    navierStokesProblem.addDirichletBC (dolfin::DirichletBC (*V[0], noSlipDirichletBC, noSlipBoundary, "topological", false));
    navierStokesProblem.addDirichletBC (dolfin::DirichletBC (*(*V[0])[1], symmetryDirichletBC, gammaSD));
*/

// NB: from FEniCS 1.5.0 C++ documentation for dolfin::DirichletBC
// The boundary indicators may be specified in a number of different ways.
// The simplest approach is to specify the boundary by a SubDomain object, using the inside() function to specify 
// on which facets the boundary conditions should be applied. The boundary facets will then be searched for and 
// marked only on the first call to apply.
// This means that the mesh could be moved after the first apply and the boundary markers would still remain intact.
// Alternatively, the boundary may be specified by a MeshFunction over facets labeling all mesh facets together
// with a number that specifies which facets should be included in the boundary.
// The third option is to attach the boundary information to the mesh.
// This is handled automatically when exporting a mesh from for example VMTK.
    myNavierstokesTimeCurv::LeftBoundary wallLeft;
    myNavierstokesTimeCurv::RightBoundary wallRight;
    Ivan::TopBoundary freeSurface;
 //   myNavierstokesTimeCurv::BottomBoundary inflowBoundary;
    dcp::Subdomain inflowBoundary ((myNavierstokesTimeCurv::BottomBoundaryEvaluator ()));
	
//sistemare ++
myNavierstokesTimeCurv::TriplePoints triplePoints;
dolfin::BoundaryMesh bdMesh (* navierStokesProblem.mesh(), "exterior");
//dolfin::plot(bdMesh.entity_map(2)); dolfin::interactive();
auto bla1 (bdMesh.entity_map(1));
auto bla0 (bdMesh.entity_map(0));
std::cerr << bla1.str(true) << "\n ciao\n" << bla0.str(true) << std::endl;
dolfin::plot(bla1); dolfin::interactive();
dolfin::plot(bla0); dolfin::interactive(); //return 0;

    //myNavierstokesTimeCurv::InflowDirichletBC inflowDirichletBC;
    dcp::TimeDependentExpression inflowDirichletBC(2,(myNavierstokesTimeCurv::InflowDirichletBCEvaluator ()));

 //   navierStokesProblem.addDirichletBC (dolfin::DirichletBC (*V[0], inflowDirichletBC, inflowBoundary));
    navierStokesProblem.addTimeDependentDirichletBC (inflowDirichletBC, inflowBoundary, 0);
//    navierStokesProblem.addDirichletBC (dolfin::DirichletBC (*V[0], noSlipDirichletBC, noSlipBoundaryTop, "topological", false));
/*    navierStokesProblem.addDirichletBC (dolfin::DirichletBC (*V[0], noSlipDirichletBC, wallLeft));
    navierStokesProblem.addDirichletBC (dolfin::DirichletBC (*V[0], noSlipDirichletBC, wallRight));
*/    navierStokesProblem.addDirichletBC (dolfin::DirichletBC (*(*V[0])[0], symmetryDirichletBC, wallLeft));
    navierStokesProblem.addDirichletBC (dolfin::DirichletBC (*(*V[0])[0], symmetryDirichletBC, wallRight));
//    navierStokesProblem.addDirichletBC (dolfin::DirichletBC (*(*V[0])[1], symmetryDirichletBC, gammaSD));

//sistemare ++
/*	dolfin::Constant tpDirichletBC (0,1);
	navierStokesProblem.addDirichletBC (dolfin::DirichletBC (*V[0], tpDirichletBC, triplePoints, "pointwise"));*/
	//CONTINUARE dolfin::VertexFunction ...
    
    // define neumann boundary conditions 
    dolfin::info ("Define the problem's Neumann boundary conditions...");
//    poisson::NeumannBoundary neumannBoundary;
    dolfin::FacetFunction<std::size_t> meshFacets (* navierStokesProblem.mesh());
    meshFacets.set_all (0);
    freeSurface.mark (meshFacets, 1);
inflowBoundary.mark (meshFacets, 2);
//    gammaSD.mark(meshFacets, 2);
    navierStokesProblem.setIntegrationSubdomain ("residual_form",
                                            dolfin::reference_to_no_delete_pointer (meshFacets),
                                            dcp::SubdomainType::BOUNDARY_FACETS);
meshFacets.set_value(10,3);
dolfin::plot(meshFacets);dolfin::interactive();
dolfin::VertexFunction<std::size_t> meshPoints (* navierStokesProblem.mesh());
meshPoints.set_all (0);
//meshPoints.set_value(10,3);
//triplePoints.mark(meshPoints,8);
dolfin::plot(meshPoints);dolfin::interactive();
Ivan::TriplePointLeft triplePointLeft;
triplePointLeft.mark(meshFacets,4);
Ivan::TriplePointRight triplePointRight;
triplePointRight.mark(meshFacets,5);
dolfin::plot(meshFacets, "meshFacets finale");dolfin::interactive();
navierStokesProblem.setIntegrationSubdomain ("residual_form",
                                            dolfin::reference_to_no_delete_pointer (meshPoints),
                                            dcp::SubdomainType::VERTICES);
navierStokesProblem.setIntegrationSubdomain ("residual_form",
                                            dolfin::reference_to_no_delete_pointer (meshFacets),
                                            dcp::SubdomainType::BOUNDARY_FACETS);
//return 0;

    // problem settings
    std::cout << "Set the problem's coefficients..." << std::endl;
    navierStokesProblem.setCoefficient ("residual_form", dolfin::reference_to_no_delete_pointer (nu), "nu");
    navierStokesProblem.setCoefficient ("jacobian_form", dolfin::reference_to_no_delete_pointer (nu), "nu");
//    navierStokesProblem.setCoefficient ("residual_form", dolfin::reference_to_no_delete_pointer (p_out), "p_out");
//    navierStokesProblem.setCoefficient ("residual_form", dolfin::reference_to_no_delete_pointer (p_left), "p_left");
    navierStokesProblem.setCoefficient ("residual_form", dolfin::reference_to_no_delete_pointer (gamma), "gamma");
myNavierstokesTimeCurv::CoefficientSpace_deltaDirac deltaFS (* navierStokesProblem.mesh());
Ivan::DeltaDirac deltaExpr;
dolfin::Function deltaFun (deltaFS);
deltaFun = deltaExpr;
    navierStokesProblem.setCoefficient ("residual_form", dolfin::reference_to_no_delete_pointer (deltaFun), "deltaDirac");
    
    navierStokesProblem.parameters ["time_stepping_solution_component"] = 0;
    
    navierStokesProblem.setInitialSolution (dolfin::Constant (0,0,0));

    // solve problem
    std::cout << "Solve the problem..." << std::endl;
    navierStokesProblem.solve ();

    // saves and plots
    dolfin::File solutionFile ("myNsCurv.pvd");
    const std::vector<std::pair<double, dolfin::Function> >& solutionsVector (navierStokesProblem.solutions()); 
    double t=t0;
    for (auto it = solutionsVector.begin(); it != solutionsVector.end(); it++, t+=dt)
 //       solutionFile << std::pair<const dolfin::Function*, double>(&(*it),t);
        solutionFile << std::pair<const dolfin::Function*, double>(&(it->second),it->first);
    dolfin::plot (* navierStokesProblem.solution().function_space()->mesh(), "Final mesh");
    dolfin::interactive ();
    dolfin::plot (* solutionsVector[5].second.function_space()->mesh(), "mesh intermedia e' uguale a finale");
    
    dolfin::interactive ();
    
    return 0;
}
