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
#include <dcp/differential_problems/differential_problems.h>
#include "poisson.h"
#include "auxForms.h"

namespace poisson
{
    class DirichletBoundary : public dolfin::SubDomain
    {
        bool inside (const dolfin::Array<double>& x, bool on_boundary) const
        {
            return dolfin::near (x[1], 0) && on_boundary;
        } 
    };
    
    class NeumannBoundary : public dolfin::SubDomain
    {
        bool inside (const dolfin::Array<double>& x, bool on_boundary) const
        {
            return false;
        } 
    };
    
    class PenaltyBoundary : public dolfin::SubDomain
    {
        bool inside (const dolfin::Array<double>& x, bool on_boundary) const
        {
            return (dolfin::near (x[1], 1) && on_boundary);
        }
    };
    
    class DatumExpr : public dolfin::Expression
    {
        void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x) const
        {
            values[0] = sin(9.42*x[0]);
        }
    };

    class ExprVec : public dolfin::Expression
    {
        void eval(dolfin::Array<double>& values, const dolfin::Array<double>& x) const
        {
            values[0] = x[0]*x[1]*x[1];
            values[1] = exp(0.5*x[0])*x[1];
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
    // create mesh and finite element space 
    dolfin::info ("Create mesh and finite element space...");
    dolfin::UnitSquareMesh mesh (20, 20);
    poisson::FunctionSpace V (mesh);
    
    // define problem
    dolfin::info ("Define the problem...");
    dcp::LinearProblem <poisson::BilinearForm, poisson::LinearForm> poissonProblem (dolfin::reference_to_no_delete_pointer (V));

    // define coefficients
    dolfin::info ("Define the problem's coefficients...");
    dolfin::Constant k (1.0);
    dolfin::Constant f (1.0);
    dolfin::Constant g (1.0);
    dolfin::Constant h (0.0);
    poisson::LinearForm::CoefficientSpace_datum datumSpace (mesh);
    dolfin::Function datum (datumSpace);
    datum = (poisson::DatumExpr());

    // define dirichlet boundary conditions 
    dolfin::info ("Define the problem's Dirichlet boundary conditions...");
    poisson::DirichletBoundary dirichletBoundary;
    poissonProblem.addDirichletBC (h, dirichletBoundary);
    
    // define neumann boundary conditions 
    dolfin::info ("Define the problem's Neumann boundary conditions...");
    poisson::NeumannBoundary neumannBoundary;
    dolfin::FacetFunction<std::size_t> meshFacets (mesh);
    meshFacets.set_all (0);
    neumannBoundary.mark (meshFacets, 1);
    poisson::PenaltyBoundary penaltyBoundary;
    penaltyBoundary.mark (meshFacets, 2);
    poissonProblem.setIntegrationSubdomain ("linear_form", 
                                             dolfin::reference_to_no_delete_pointer (meshFacets), 
                                             dcp::SubdomainType::BOUNDARY_FACETS);
    poissonProblem.setIntegrationSubdomain ("bilinear_form", 
                                             dolfin::reference_to_no_delete_pointer (meshFacets), 
                                             dcp::SubdomainType::BOUNDARY_FACETS);
    
    // set problem coefficients
    dolfin::info ("Set the problem's coefficients...");
    poissonProblem.setCoefficient ("bilinear_form", dolfin::reference_to_no_delete_pointer (k), "k");
    poissonProblem.setCoefficient ("linear_form", dolfin::reference_to_no_delete_pointer (f), "f");
    poissonProblem.setCoefficient ("linear_form", dolfin::reference_to_no_delete_pointer (g), "g");
    poissonProblem.setCoefficient ("linear_form", dolfin::reference_to_no_delete_pointer (datum), "datum");

    // solve problem
    dolfin::info ("Solve the problem...");
    poissonProblem.solve ();
    
    // plots
    dolfin::plot (mesh);
    poissonProblem.plotSolution ();

    // forms check
    auxForms::Form_diffP1::TestSpace diffV (mesh);
    auxForms::Form_diffP2::TestSpace diffV2 (mesh);
    auxForms::Form_diffP4::TestSpace diffV4 (mesh);
    auxForms::Form_diffP6::TestSpace diffV6 (mesh);
    dolfin::Function fun4 (diffV4);
    dolfin::Function fun6 (diffV6);
    auxForms::Form_diffP1 diffForm (diffV);
    auxForms::Form_diffP2 diffForm2 (diffV2);
    auxForms::Form_diffP4 diffForm4 (diffV4);
    auxForms::Form_diffP6 diffForm6 (diffV6);
    dolfin::Vector b (* poissonProblem.solution().vector());
    dolfin::Vector b2 (* datum.vector());
    dolfin::Vector b4 (* fun4.vector());
    dolfin::Vector b6 (* fun6.vector());

    auxForms::Form_diffP1::CoefficientSpace_coeff coeffSpace (mesh);
    dolfin::Function coeff (coeffSpace);
    std::srand(10);
    for (std::size_t i (0); i<coeff.vector()->size(); ++i)
    {
        coeff.vector()->setitem (i, 0.2/RAND_MAX*std::rand());
        coeff.vector()->apply("insert");
    }

    diffForm.set_coefficient ("u", dolfin::reference_to_no_delete_pointer (poissonProblem.solution()));
    diffForm.set_coefficient ("datum", dolfin::reference_to_no_delete_pointer (datum));
    diffForm.set_coefficient ("coeff", dolfin::reference_to_no_delete_pointer (coeff));
    diffForm2.set_coefficient ("u", dolfin::reference_to_no_delete_pointer (poissonProblem.solution()));
    diffForm2.set_coefficient ("datum", dolfin::reference_to_no_delete_pointer (datum));
    diffForm2.set_coefficient ("coeff", dolfin::reference_to_no_delete_pointer (coeff));
    diffForm4.set_coefficient ("u", dolfin::reference_to_no_delete_pointer (poissonProblem.solution()));
    diffForm4.set_coefficient ("datum", dolfin::reference_to_no_delete_pointer (datum));
    diffForm4.set_coefficient ("coeff", dolfin::reference_to_no_delete_pointer (coeff));
    diffForm6.set_coefficient ("u", dolfin::reference_to_no_delete_pointer (poissonProblem.solution()));
    diffForm6.set_coefficient ("datum", dolfin::reference_to_no_delete_pointer (datum));
    diffForm6.set_coefficient ("coeff", dolfin::reference_to_no_delete_pointer (coeff));
    diffForm.set_exterior_facet_domains (dolfin::reference_to_no_delete_pointer (meshFacets));
    diffForm2.set_exterior_facet_domains (dolfin::reference_to_no_delete_pointer (meshFacets));
    diffForm4.set_exterior_facet_domains (dolfin::reference_to_no_delete_pointer (meshFacets));
    diffForm6.set_exterior_facet_domains (dolfin::reference_to_no_delete_pointer (meshFacets));
    dolfin::Assembler assembler;
    b.zero();
    b2.zero();
    b4.zero();
    b6.zero();
    std::cerr << "                      L1,  \t\tL2,\t\tL_infty" << std::endl;
    std::cerr << "random coefficient (P2)" << std::endl;
    assembler.assemble (b, diffForm);
    std::cerr << " error (test P1) = " << b.sum() << ",  " << b.inner(b) << ",  " << b.max() << std::endl;
    assembler.assemble (b2, diffForm2);
    std::cerr << " error (test P2) = " << b2.sum() << ",  " << b2.inner(b2) << ",  " << b2.max() << std::endl;
    assembler.assemble (b4, diffForm4);
    std::cerr << " error (test P4) = " << b4.sum() << ",  " << b4.inner(b4) << ",  " << b4.max() << std::endl;
    assembler.assemble (b6, diffForm6);
    std::cerr << " error (test P6) = " << b6.sum() << ",  " << b6.inner(b6) << ",  " << b6.max() << std::endl;

    coeff.vector()->zero();
    (* coeff.vector()) += 1.0;
    diffForm.set_coefficient ("coeff", dolfin::reference_to_no_delete_pointer (coeff));
    diffForm2.set_coefficient ("coeff", dolfin::reference_to_no_delete_pointer (coeff));
    diffForm4.set_coefficient ("coeff", dolfin::reference_to_no_delete_pointer (coeff));
    diffForm6.set_coefficient ("coeff", dolfin::reference_to_no_delete_pointer (coeff));
    assembler.assemble (b, diffForm);
    std::cerr << "unit coefficient (P2)" << std::endl;
    std::cerr << " error (test P1) = " << b.sum() << ",  " << b.inner(b) << ",  " << b.max() << std::endl;
    assembler.assemble (b2, diffForm2);
    std::cerr << " error (test P2) = " << b2.sum() << ",  " << b2.inner(b2) << ",  " << b2.max() << std::endl;
    assembler.assemble (b4, diffForm4);
    std::cerr << " error (test P4) = " << b4.sum() << ",  " << b4.inner(b4) << ",  " << b4.max() << std::endl;
    assembler.assemble (b6, diffForm6);
    std::cerr << " error (test P6) = " << b6.sum() << ",  " << b6.inner(b6) << ",  " << b6.max() << std::endl;
    
    dolfin::interactive ();
    
    dolfin::Constant unoVec (0.03,0.04);
    auxForms::Form_kinEnNO::CoefficientSpace_uVec uVecSpace (mesh);
    dolfin::Function uVec (uVecSpace);
    uVec = unoVec;

    auxForms::Form_kinEnNO::TestSpace kinEnNOSpace (mesh);
    auxForms::Form_kinEn::TestSpace kinEnSpace (mesh);
    auxForms::Form_kinEn2::TestSpace kinEn2Space (mesh);
    auxForms::Form_kinEnNO kinEnNO (kinEnNOSpace);
    auxForms::Form_kinEn kinEn (kinEnSpace);
    auxForms::Form_kinEn2 kinEn2 (kinEn2Space);

    kinEnNO.set_coefficient ("uVec", dolfin::reference_to_no_delete_pointer (uVec));
    kinEn.set_coefficient ("uVec", dolfin::reference_to_no_delete_pointer (uVec));
    kinEn2.set_coefficient ("uVec", dolfin::reference_to_no_delete_pointer (uVec));
    dolfin::Vector eNO (MPI_COMM_WORLD, kinEnNOSpace.dim());
    dolfin::Vector e (MPI_COMM_WORLD, kinEnSpace.dim());
    dolfin::Vector e2 (MPI_COMM_WORLD, kinEn2Space.dim());
    assembler.assemble (eNO, kinEnNO);
    assembler.assemble (e, kinEn);
    assembler.assemble (e2, kinEn2);
    std::cerr << "sums " << eNO.sum() << ", " << e.sum() << ", " << e2.sum() << std::endl;

    uVec = (poisson::ExprVec());
    kinEnNO.set_coefficient ("uVec", dolfin::reference_to_no_delete_pointer (uVec));
    kinEn.set_coefficient ("uVec", dolfin::reference_to_no_delete_pointer (uVec));
    kinEn2.set_coefficient ("uVec", dolfin::reference_to_no_delete_pointer (uVec));
    assembler.assemble (eNO, kinEnNO);
    assembler.assemble (e, kinEn);
    assembler.assemble (e2, kinEn2);
    std::cerr << "sums " << eNO.sum() << ", " << e.sum() << ", " << e2.sum() << std::endl;

    return 0;
}
