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
#include <mshr.h>
#include "primal.h"
#include "lift_drag.h"
#include <differential_problems/NonlinearDifferentialProblem.h>

namespace navierstokes
{
    class InflowBoundary : public dolfin::SubDomain
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
}

int main (int argc, char* argv[])
{
    // define parameters and their default values
    dolfin::Parameters parameters ("main_parameters");
    parameters.add ("output_file_name", "target_solution");
    parameters.add ("human_readable_print", false);
   
    // read parameters from command line and overwrite default values
    parameters.parse (argc, argv);
    
    // create mesh and finite element space
    mshr::Rectangle domain (dolfin::Point (0.0, 0.0), dolfin::Point (10.0, 7.0));
    mshr::Rectangle rectangle (dolfin::Point (1.75, 3.0), dolfin::Point (2.0, 4.0));
    mshr::Circle obstacle (dolfin::Point (3.5, 3.5), 0.5);
    dolfin::Mesh mesh;
    mshr::generate (mesh, *(domain - rectangle - obstacle), 40);
    primal::FunctionSpace V (mesh);
    
    
    // define constant
    dolfin::Constant nu (1e-1);
    dolfin::Constant inflowDirichletBC (1.0, 0.0);
    dolfin::Constant symmetryDirichletBC (0.0);
    dolfin::Constant noSlipDirichletBC (0.0, 0.0);

    // define boundary conditions subdomains
    navierstokes::InflowBoundary inflowBoundary;
    navierstokes::GammaSD gammaSD;
    navierstokes::NoSlipBoundary noSlipBoundary;

    // define problem
    dcp::NonlinearDifferentialProblem <primal::ResidualForm, primal::JacobianForm> 
        navierStokesProblem (dolfin::reference_to_no_delete_pointer (mesh), 
                       dolfin::reference_to_no_delete_pointer (V),
                       "trial");

    // problem settings
    navierStokesProblem.setCoefficient ("residual_form", dolfin::reference_to_no_delete_pointer (nu), "nu");
    navierStokesProblem.setCoefficient ("jacobian_form", dolfin::reference_to_no_delete_pointer (nu), "nu");

    navierStokesProblem.addDirichletBC (dolfin::DirichletBC (*V[0], inflowDirichletBC, inflowBoundary));
    navierStokesProblem.addDirichletBC (dolfin::DirichletBC (*V[0], noSlipDirichletBC, noSlipBoundary, "topological", false));
    navierStokesProblem.addDirichletBC (dolfin::DirichletBC (*(*V[0])[1], symmetryDirichletBC, gammaSD));

    
    // solve problem
    navierStokesProblem.solve ();

    // plots
    dolfin::plot (mesh, "Mesh");
    dolfin::plot (navierStokesProblem.solution ()[0], "Velocity");
    dolfin::plot (navierStokesProblem.solution ()[1], "Pressure");
    dolfin::interactive ();
    
    // print to file
    dolfin::cout << "Printing to hdf5 format..." << dolfin::endl;
    dolfin::HDF5File uFile (MPI_COMM_WORLD, static_cast<std::string> (parameters ["output_file_name"]) + ".hdf5", "w");
    uFile.write (navierStokesProblem.solution (), "solution");
    
    if (static_cast<bool> (parameters ["human_readable_print"]) == true)
    {    
        dolfin::cout << "Printing to human readable format..." << dolfin::endl;
        dolfin::File uPlotFile (static_cast<std::string> (parameters ["output_file_name"]) + "_velocity" + ".pvd");
        dolfin::Function u (navierStokesProblem.solution ()[0]);
        uPlotFile << u;
        
        dolfin::File pPlotFile (static_cast<std::string> (parameters ["output_file_name"]) + "_pressure" + ".pvd");
        dolfin::Function p (navierStokesProblem.solution ()[1]);
        pPlotFile << p;
    }
    
    // ============================================================================== //
    // =============================== LIFT AND DRAG ================================ //
    // ============================================================================== //
    dolfin::FacetFunction<std::size_t> meshFacets (mesh);
    meshFacets.set_all (0);
    
    navierstokes::Circle circle;
    circle.mark (meshFacets, 3);
    
    // ------
    // lift
    // ------
    lift_drag::Form_Lift liftComputer (mesh);
    
    // functional settings
    liftComputer.set_coefficient ("p", dolfin::reference_to_no_delete_pointer (navierStokesProblem.solution ()[1]));
    liftComputer.set_exterior_facet_domains (dolfin::reference_to_no_delete_pointer (meshFacets));
    
    
    // ------
    // drag
    // ------
    lift_drag::Form_Drag dragComputer (mesh);
        
    // functional settings
    dragComputer.set_coefficient ("p", dolfin::reference_to_no_delete_pointer (navierStokesProblem.solution ()[1]));
    dragComputer.set_exterior_facet_domains (dolfin::reference_to_no_delete_pointer (meshFacets));
    
   
    // compute lift and drag
    std::ofstream liftDragOutputStream ("target_lift_drag.txt");
    liftDragOutputStream << "Lift = " << dolfin::assemble (liftComputer) << std::endl;
    liftDragOutputStream << "Drag = " << dolfin::assemble (dragComputer) << std::endl;
    liftDragOutputStream.close ();
    
    return 0;
}
