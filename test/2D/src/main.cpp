#include <iostream>
#include <dolfin.h>
#include "primal.h"
#include "adjoint.h"
#include <DifferentialProblem/NonlinearDifferentialProblem.hpp>
#include <DifferentialProblem/LinearDifferentialProblem.hpp>
#include <DifferentialProblem/CompositeDifferentialProblem.hpp>

namespace primal
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
            return (x[1] < (0 + DOLFIN_EPS) || x[1] > (2 - DOLFIN_EPS)) && on_boundary;
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
                double dx = x[0] - center[0];
                double dy = x[1] - center[1];
                double r = sqrt (dx * dx + dy * dy);
                
                return r < (radius + 1e-3) && on_boundary;
            }

        private:
            dolfin::Array<double> center;
            double radius;
    }; 
}

namespace adjoint
{
    class DirichletBoundary : public dolfin::SubDomain
    {
        public:
            DirichletBoundary () : 
                center (2),
                radius (0.5)
            {
                center[0] = 2.5;
                center[1] = 1;
            }

            bool inside (const dolfin::Array<double>& x, bool on_boundary) const
            {
                bool onLeftSide  = (x[0] < (0 + DOLFIN_EPS)) && on_boundary;   
                bool onUpperSide = (x[1] > (2 - DOLFIN_EPS)) && on_boundary;   
                bool onLowerSide = (x[1] < (0 + DOLFIN_EPS)) && on_boundary;   

                double dx = x[0] - center[0];
                double dy = x[1] - center[1];
                double r = sqrt (dx * dx + dy * dy);
                bool onCircleBoundary = r < (radius + 1e-3) && on_boundary;

                return onLeftSide || onUpperSide || onLowerSide || onCircleBoundary;
            }
        
        private:
            dolfin::Array<double> center;
            double radius;
    };

    class RobinBoundary : public dolfin::SubDomain
    {
        bool inside (const dolfin::Array<double>& x, bool on_boundary) const
        {
            return x[0] > (5 - DOLFIN_EPS) && on_boundary;
        }
    };
    
    class ExternalLoadDomain : public dolfin::SubDomain
    {
        bool inside (const dolfin::Array<double>& x, bool on_boundary) const
        {
            return x[0] >= 1.75 && x[0] <= 3.25 && x[1] >= 0 && x[1] <= 1.75; 
        }
    };
}


int main (int argc, char* argv[])
{
    // log level
//    dolfin::set_log_level (dolfin::DBG);

    // define parameters and their default values
    dolfin::Parameters parameters ("main_parameters");
    parameters.add ("complete_mesh_file_name", "../src/complete_mesh/complete_mesh.xml");
    parameters.add ("partial_mesh_file_name", "../src/partial_mesh/partial_mesh.xml");
    parameters.add ("target_u_file_name", "../src/u_target");
    parameters.add ("target_p_file_name", "../src/p_target");
    
    // define parameter to allow projection of solution from complete mesh to partial mesh.
    // In any case, the values that will be extrapolated are not important because the functional will
    // integrate the difference of the two only on a subdomain, which does not include the points whre
    // extrapolation takes place
    dolfin::parameters ["allow_extrapolation"] = true;
    
    // read parameters from command line and overwrite default values
    parameters.parse (argc, argv);
    
    // create mesh and finite element space to read target solution
    dolfin::Mesh completeMesh (parameters ["complete_mesh_file_name"]);
    primal::FunctionSpace completeV (completeMesh);
    
    // read target solutions from file
    dolfin::Function targetU (completeV[0] -> collapse ());
    dolfin::Function targetP (completeV[1] -> collapse ());
    
    dolfin::HDF5File targetUFile (static_cast<std::string> (parameters ["target_u_file_name"]) + ".hdf5", "r");
    targetUFile.read (targetU, "u");
    dolfin::HDF5File targetPFile (static_cast<std::string> (parameters ["target_p_file_name"]) + ".hdf5", "r");
    targetPFile.read (targetP, "p");
    
    
    // create mesh and finite element space
    dolfin::Mesh mesh (parameters ["partial_mesh_file_name"]);
    primal::FunctionSpace V (mesh);
    
    dolfin::plot (mesh);
    
    
    // interpolate target functions from completeMesh to actual mesh
    dolfin::Function interpolatedTargetU (V[0] -> collapse ());
    dolfin::Function interpolatedTargetP (V[1] -> collapse ());
    
    interpolatedTargetU.interpolate (targetU);
    interpolatedTargetP.interpolate (targetP);

    
    // define problems
    controlproblem::NonlinearDifferentialProblem <primal::ResidualForm, primal::JacobianForm> 
        primalProblem2 (dolfin::reference_to_no_delete_pointer (mesh), 
                        dolfin::reference_to_no_delete_pointer (V),
                        "trial");

    controlproblem::LinearDifferentialProblem <adjoint::BilinearForm, adjoint::LinearForm> 
        adjointProblem2 (dolfin::reference_to_no_delete_pointer (mesh),
                         dolfin::reference_to_no_delete_pointer (V));

    // create composite differential problem
    controlproblem::CompositeDifferentialProblem problems;
    problems.addProblem ("primal", primalProblem2);
    problems.addProblem ("adjoint", adjointProblem2);
    
    // define constants
    dolfin::Constant nu (1e-5);
    dolfin::Constant primal_inflowDirichletBC (1.0, 0.0);
    dolfin::Constant primal_symmetryDirichletBC (0.0);
    dolfin::Constant primal_noSlipCondition (0.0, 0.0);
    dolfin::Constant adjoint_dirichletBC (0.0, 0.0);

    // define boundary conditions subdomains
    primal::InflowBoundary primal_inflowBoundary;
    primal::GammaSD primal_gammaSD;
    primal::NoSlipBoundary primal_noSlipBoundary;
    adjoint::DirichletBoundary adjoint_dirichletBoundary;
    adjoint::RobinBoundary adjoint_robinBoundary;
    adjoint::ExternalLoadDomain adjoint_externalLoadDomain;
    
    // define intergration subdomains
    dolfin::FacetFunction<std::size_t> adjoint_meshFacets (mesh);
    adjoint_meshFacets.set_all (1);
    adjoint_robinBoundary.mark (adjoint_meshFacets, 0);
    
    dolfin::CellFunction<std::size_t> adjoint_meshCells (mesh);
    adjoint_meshCells.set_all (1);
    adjoint_externalLoadDomain.mark (adjoint_meshCells, 0);

    // primal problem settings
    problems["primal"].setCoefficient ("residual_form", dolfin::reference_to_no_delete_pointer (nu), "nu");
    problems["primal"].setCoefficient ("jacobian_form", dolfin::reference_to_no_delete_pointer (nu), "nu");

    problems["primal"].addDirichletBC (dolfin::DirichletBC (*V[0], primal_inflowDirichletBC, primal_inflowBoundary));
    problems["primal"].addDirichletBC (dolfin::DirichletBC (*V[0], primal_noSlipCondition, primal_noSlipBoundary));
    problems["primal"].addDirichletBC (dolfin::DirichletBC (*(*V[0])[1], primal_symmetryDirichletBC, primal_gammaSD));


    // adjoint problem settings
    problems["adjoint"].setCoefficient ("bilinear_form", dolfin::reference_to_no_delete_pointer (nu), "nu");
    problems["adjoint"].setCoefficient ("linear_form", dolfin::reference_to_no_delete_pointer (interpolatedTargetU), "U");
    problems["adjoint"].setCoefficient ("linear_form", dolfin::reference_to_no_delete_pointer (interpolatedTargetP), "P");

    problems["adjoint"].addDirichletBC (dolfin::DirichletBC (*V[0], adjoint_dirichletBC, adjoint_dirichletBoundary));

    problems["adjoint"].setIntegrationSubdomains ("bilinear_form", 
                                                  dolfin::reference_to_no_delete_pointer (adjoint_meshFacets),
                                                  controlproblem::SubdomainType::BOUNDARY_FACETS);
    
    problems["adjoint"].setIntegrationSubdomains ("linear_form", 
                                                  dolfin::reference_to_no_delete_pointer (adjoint_meshCells),
                                                  controlproblem::SubdomainType::INTERNAL_CELLS);

    problems.addLink ("adjoint", "u", "bilinear_form", "primal", 0);
    problems.addLink ("adjoint", "u", "linear_form", "primal", 0);
    problems.addLink ("adjoint", "p", "linear_form", "primal", 1);
    
    problems.print ();
    
    
    // solve problems
    problems.solve ();

    // plots
    dolfin::plot (problems.solution ("primal")[0]);
    dolfin::plot (problems.solution ("primal")[1]);
    dolfin::plot (problems.solution ("adjoint")[0]);
    dolfin::plot (problems.solution ("adjoint")[1]);

    dolfin::interactive ();

    return 0;
}
