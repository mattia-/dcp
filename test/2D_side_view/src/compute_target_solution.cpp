#include <iostream>
#include <string>
#include <dolfin.h>
#include "primal.h"
#include <DifferentialProblem/NonlinearDifferentialProblem.hpp>

namespace navierstokes
{
    class InflowBoundary : public dolfin::SubDomain
    {
        bool inside (const dolfin::Array<double>& x, bool on_boundary) const
        {
            return x[0] < (0 + DOLFIN_EPS) && on_boundary;
        }
    };

    class GammaSD : public dolfin::SubDomain
    {
        bool inside (const dolfin::Array<double>& x, bool on_boundary) const
        {
            return (x[1] < (0 + DOLFIN_EPS) || x[1] > (2.5 - DOLFIN_EPS)) && on_boundary;
        }
    };
    
    class NoSlipBoundary : public dolfin::SubDomain
    { 
        public:
            NoSlipBoundary () : 
                center (2),
                radius (0.5)
            {
                center[0] = 2.5;
                center[1] = 1;
            }
            
            bool inside (const dolfin::Array<double>& x, bool on_boundary) const
            {
                bool onWall1  = (x[0] > (1.25 - DOLFIN_EPS) && x[0] < (1.25 + DOLFIN_EPS))
                    && x[1] < (1.5 + DOLFIN_EPS)
                    && on_boundary;
                bool onWall2  = (x[0] > (1.25 - DOLFIN_EPS) && x[0] < (1.5 + DOLFIN_EPS))
                    && (x[1] < (1.5 + DOLFIN_EPS) && x[1] > (1.5 - DOLFIN_EPS))
                    && on_boundary;
                bool onWall3  = (x[0] > (1.5 - DOLFIN_EPS) && x[0] < (1.5 + DOLFIN_EPS))
                    && x[1] < (1.5 + DOLFIN_EPS)
                    && on_boundary;
                
                double dx = x[0] - center[0];
                double dy = x[1] - center[1];
                double r = sqrt (dx * dx + dy * dy);
                bool onCircleBoundary = r < (radius + 1e-3) && on_boundary;

                return onWall1 || onWall2 || onWall3 || onCircleBoundary;
            }

        private:
            dolfin::Array<double> center;
            double radius;
    }; 
}

int main (int argc, char* argv[])
{
    // define parameters and their default values
    dolfin::Parameters parameters ("main_parameters");
    parameters.add ("mesh_file_name", "../src/complete_mesh/complete_mesh.xml");
    parameters.add ("u_output_file_name", "../src/u_target");
    parameters.add ("p_output_file_name", "../src/p_target");
    parameters.add ("human_readable_print", false);
    
    // read parameters from command line and overwrite default values
    parameters.parse (argc, argv);
    
    // create mesh and finite element space
    dolfin::Mesh mesh (parameters ["mesh_file_name"]);
    primal::FunctionSpace V (mesh);
    
    dolfin::plot (mesh);

    // define constant
    dolfin::Constant nu (1e-2);
    dolfin::Constant inflowDirichletBC (1.0, 0.0);
    dolfin::Constant symmetryDirichletBC (0.0);
    dolfin::Constant noSlipDirichletBC (0.0, 0.0);

    // define boundary conditions subdomains
    navierstokes::InflowBoundary inflowBoundary;
    navierstokes::GammaSD gammaSD;
    navierstokes::NoSlipBoundary noSlipBoundary;

    // define problem
    controlproblem::NonlinearDifferentialProblem <primal::ResidualForm, primal::JacobianForm> 
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
    dolfin::plot (navierStokesProblem.solution ()[0]);
    dolfin::plot (navierStokesProblem.solution ()[1]);

    dolfin::interactive ();
    
    
    // print to file
    dolfin::cout << "Printing to hdf5 format..." << dolfin::endl;
    dolfin::HDF5File uFile (static_cast<std::string> (parameters ["u_output_file_name"]) + ".hdf5", "w");
    uFile.write (navierStokesProblem.solution ()[0], "u");
    
    dolfin::HDF5File pFile (static_cast<std::string> (parameters ["p_output_file_name"]) + ".hdf5", "w");
    pFile.write (navierStokesProblem.solution ()[1], "p");

    if (static_cast<bool> (parameters ["human_readable_print"]) == true)
    {    
        dolfin::cout << "Printing to human readable format..." << dolfin::endl;
        dolfin::File uPlotFile (static_cast<std::string> (parameters ["u_output_file_name"]) + ".pvd");
        dolfin::Function u (navierStokesProblem.solution ()[0]);
        uPlotFile << u;
        
        dolfin::File pPlotFile (static_cast<std::string> (parameters ["p_output_file_name"]) + ".pvd");
        dolfin::Function p (navierStokesProblem.solution ()[1]);
        pPlotFile << p;
    }
    
    return 0;
}
