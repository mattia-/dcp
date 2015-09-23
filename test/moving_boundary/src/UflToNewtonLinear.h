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

#ifndef IVAN_UFLTONEWTONLINEAR_H_INCLUDE_GUARD
#define IVAN_UFLTONEWTONLINEAR_H_INCLUDE_GUARD

#include <dolfin.h>
#include <dcp/differential_problems/SubdomainType.h>
#include <dcp/differential_problems/LinearProblem.h>
#include "utilities.h"

namespace Ivan
{

template <class T_FunctionSpace, class T_AdditionalFunctionSpace, class T_BilinearForm, class T_LinearForm, class T_AdditionalForm, class T_LinearSolverFactory = dcp::LinearSolverFactory>
class UflToNewtonLinear : public dcp::LinearProblem <T_BilinearForm, T_LinearForm, T_LinearSolverFactory>
{

  public:

    // Constructor
/*    UflToNewtonLinear (const std::shared_ptr<const dolfin::VertexFunction<std::size_t> > meshFunction,
                       const std::set<std::size_t> labelSet,
                       std::shared_ptr<dolfin::FunctionSpace> functionSpace);*/
    UflToNewtonLinear (std::shared_ptr<dolfin::FunctionSpace> functionSpace, const std::vector<dolfin::la_index> & additionalFormDofs);

    // Return solution function (non const version)
    // TODO sarebbe meglio tenere solo quella const di dcp::LinearProblem;
    //      per il momento questa serve per poter passare *solution_.vector() a dolfin::LinearSolver
    dolfin::Function& solution ();

    // Overriding dcp::LinearProblem::setCoefficient [1]
    virtual void setCoefficient (const std::string& coefficientType, 
                                 const std::shared_ptr<const dolfin::GenericFunction> coefficientValue,
                                 const std::string& coefficientName = "default") override;

    // Overriding dcp::LinearProblem::setCoefficient [2]
    virtual void setCoefficient (const std::string& coefficientType,
                                 const std::shared_ptr<const dolfin::GenericFunction> coefficientValue,
                                 const std::size_t& coefficientNumber) override;

    // Overriding dcp::LinearProblem::setIntegrationSubdomain
    virtual void setIntegrationSubdomain (const std::string& formType,
                                           std::shared_ptr<const dolfin::MeshFunction<std::size_t>> meshFunction,
                                           const dcp::SubdomainType& subdomainType) override;

    // Set solver for algebraic solution
    virtual bool setLinearSolver (std::shared_ptr<dolfin::GenericLinearSolver> solver);

    //! Solve problem
    /*!
     *  This method solves the problem defined. It uses the private members' value to set the problem and then
     *  stores the solution in the private member \c solution_. See documentation of \c dcp::AbstractProblem
     *  for more details on how the protected member \c solution_ works and why it is declared as a 
     *  \c std::pair.
     *
     *  \param type the solution type requested. In this class, the only possibility is to set
     *  \c type equal to \c "default".
     */
    virtual void solve (const std::string& type = "default") override;

    protected:

      virtual void processAdditionalVector (dolfin::Vector& processedVec, const dolfin::Vector& inputVec);

/*      const std::shared_ptr<const dolfin::VertexFunction<std::size_t> > meshFunction_;
      const std::set<std::size_t> labelSet_;*/
      const std::vector<dolfin::la_index> & additionalFormDofs_;
      dolfin::FunctionSpace additionalFunctionSpace_;
      T_AdditionalForm additionalForm_;

      std::shared_ptr<dolfin::GenericLinearSolver> linear_solver_;
      bool linear_solver_set_;
  
}; //end of class definition


    // ==============================================================================================//
    // ==================================== IMPLEMENTATION ==========================================//
    // ==============================================================================================//

template <class T_FunctionSpace, class T_AdditionalFunctionSpace, class T_BilinearForm, class T_LinearForm, class T_AdditionalForm, class T_LinearSolverFactory>
UflToNewtonLinear <T_FunctionSpace, T_AdditionalFunctionSpace, T_BilinearForm, T_LinearForm, T_AdditionalForm, T_LinearSolverFactory>::
/*  UflToNewtonLinear (const std::shared_ptr<const dolfin::VertexFunction<std::size_t> > meshFunction,
                     const std::set<std::size_t> labelSet,
                     std::shared_ptr<dolfin::FunctionSpace> functionSpace) :*/
  UflToNewtonLinear (std::shared_ptr<dolfin::FunctionSpace> functionSpace, const std::vector<dolfin::la_index> & additionalFormDofs) :
    dcp::LinearProblem<T_BilinearForm,T_LinearForm,T_LinearSolverFactory> (functionSpace),
/*    meshFunction_ (meshFunction),
    labelSet_ (labelSet),*/
    additionalFormDofs_ (additionalFormDofs),
    additionalFunctionSpace_ (T_AdditionalFunctionSpace(this->functionSpace_->mesh())),
    additionalForm_ (T_AdditionalForm (this->additionalFunctionSpace_))
  {}

template <class T_FunctionSpace, class T_AdditionalFunctionSpace, class T_BilinearForm, class T_LinearForm, class T_AdditionalForm, class T_LinearSolverFactory>
dolfin::Function& UflToNewtonLinear <T_FunctionSpace, T_AdditionalFunctionSpace, T_BilinearForm, T_LinearForm, T_AdditionalForm, T_LinearSolverFactory>::
  solution ()
  {
      return this->solution_.back ().second;
  }

template <class T_FunctionSpace, class T_AdditionalFunctionSpace, class T_BilinearForm, class T_LinearForm, class T_AdditionalForm, class T_LinearSolverFactory>
void UflToNewtonLinear <T_FunctionSpace, T_AdditionalFunctionSpace, T_BilinearForm, T_LinearForm, T_AdditionalForm, T_LinearSolverFactory>::
  setCoefficient (const std::string& coefficientType, 
                  const std::shared_ptr<const dolfin::GenericFunction> coefficientValue,
                  const std::string& coefficientName)
  {
      if (coefficientType == "additional_form")
      {
          dolfin::log (dolfin::DBG, "Setting linear form coefficient \"%s\"...", coefficientName.c_str ());
          this->additionalForm_.set_coefficient (coefficientName, coefficientValue);
      }
      else
          this->dcp::LinearProblem <T_BilinearForm,T_LinearForm,T_LinearSolverFactory>::setCoefficient (coefficientType,coefficientValue,coefficientName);
  }

template <class T_FunctionSpace, class T_AdditionalFunctionSpace, class T_BilinearForm, class T_LinearForm, class T_AdditionalForm, class T_LinearSolverFactory>
void UflToNewtonLinear <T_FunctionSpace, T_AdditionalFunctionSpace, T_BilinearForm, T_LinearForm, T_AdditionalForm, T_LinearSolverFactory>::
  setCoefficient (const std::string& coefficientType,
                  const std::shared_ptr<const dolfin::GenericFunction> coefficientValue,
                  const std::size_t& coefficientNumber)
  {
      if (coefficientType == "additional_form")
      {
          dolfin::log (dolfin::DBG, "Setting linear form coefficient number %d...", coefficientNumber);
          this->additionalForm_.set_coefficient (coefficientNumber, coefficientValue);
      }
      else
          this->dcp::LinearProblem <T_BilinearForm,T_LinearForm,T_LinearSolverFactory>::setCoefficient (coefficientType,coefficientValue,coefficientNumber);
  }

template <class T_FunctionSpace, class T_AdditionalFunctionSpace, class T_BilinearForm, class T_LinearForm, class T_AdditionalForm, class T_LinearSolverFactory>
void UflToNewtonLinear <T_FunctionSpace, T_AdditionalFunctionSpace, T_BilinearForm, T_LinearForm, T_AdditionalForm, T_LinearSolverFactory>::
  setIntegrationSubdomain (const std::string& formType,
                            std::shared_ptr<const dolfin::MeshFunction<std::size_t>> meshFunction,
                            const dcp::SubdomainType& subdomainType)
  {
      if (formType == "additional_form")
      {
          if (subdomainType == dcp::SubdomainType::INTERNAL_CELLS)
          {
              dolfin::log (dolfin::DBG, "Setting linear form integration subdomain on INTERNAL_CELLS...");
              this->additionalForm_.set_cell_domains (meshFunction);
          }
          else if (subdomainType == dcp::SubdomainType::INTERNAL_FACETS)
          {
              dolfin::log (dolfin::DBG, "Setting linear form integration subdomain on INTERNAL_FACETS...");
              this->additionalForm_.set_interior_facet_domains (meshFunction);
          }
          else if (subdomainType == dcp::SubdomainType::BOUNDARY_FACETS)
          {
              dolfin::log (dolfin::DBG, "Setting linear form integration subdomain on BOUNDARY_FACETS...");
              this->additionalForm_.set_exterior_facet_domains (meshFunction);
          }
          else if (subdomainType == dcp::SubdomainType::VERTICES)
          {
              dolfin::log (dolfin::DBG, "Setting linear form integration subdomain on VERTICES...");
              this->additionalForm_.set_vertex_domains (meshFunction);
          }
          else
          {
              dolfin::warning ("unknown subdomain type requested while trying to apply mesh function to linear form"); 
          }
      }
      else
          this->dcp::LinearProblem <T_BilinearForm,T_LinearForm,T_LinearSolverFactory>::setIntegrationSubdomain (formType,meshFunction,subdomainType);
  }

template <class T_FunctionSpace, class T_AdditionalFunctionSpace, class T_BilinearForm, class T_LinearForm, class T_AdditionalForm, class T_LinearSolverFactory>
bool UflToNewtonLinear <T_FunctionSpace, T_AdditionalFunctionSpace, T_BilinearForm, T_LinearForm, T_AdditionalForm, T_LinearSolverFactory>::
  setLinearSolver (std::shared_ptr<dolfin::GenericLinearSolver> solver)
    //TODO templatizzare/altre opzioni per avere generico solutore
  {
      this->linear_solver_.reset(new dolfin::LinearSolver);
      // TODO sarebbe meglio una cosa tipo la riga seguente...
      // this->linear_solver_.reset (solver);
      linear_solver_set_ = true;
// !!! insert here possible parameters for an iterative linear solver (cfr. KrylovSolver::default_parameters())
// !!! linear_solver_->parameters["..."] = ...;
      return linear_solver_set_;
  }

template <class T_FunctionSpace, class T_AdditionalFunctionSpace, class T_BilinearForm, class T_LinearForm, class T_AdditionalForm, class T_LinearSolverFactory>
void UflToNewtonLinear <T_FunctionSpace, T_AdditionalFunctionSpace, T_BilinearForm, T_LinearForm, T_AdditionalForm, T_LinearSolverFactory>::
  solve (const std::string& type)
  {
std::cerr << "STAI USANDO L'ALGEBRICO" << std::endl;

        dolfin::Assembler assembler;

        dolfin::Matrix A;
        assembler.assemble (A, this->bilinearForm_);

        dolfin::Vector b;
        assembler.assemble (b, this->linearForm_);

        dolfin::Vector vec (MPI_COMM_WORLD,b.size());
dolfin::plot (*additionalFunctionSpace_.mesh(),"additional mesh"); //dolfin::interactive ();
//TODO sistemare: additionalFunctionSpace non subisce il mesh-moving
//     forse e' il caso di introdurre un TimeDependentFunctionSpace, cosi' tutto passa attraverso di lui...
        assembler.assemble (vec, this->additionalForm_);

        dolfin::Vector processedVec (MPI_COMM_WORLD,vec.size());
        processedVec.zero();
        processAdditionalVector (processedVec, vec);
        b += processedVec;

        for (auto it = this->dirichletBCs_.begin(); it != this->dirichletBCs_.end(); it++)
            {(it->second).apply (A,b); std::cerr<<"TimeBC "<<it->first<<std::endl;}

        linear_solver_->solve(A, * this->solution_.back().second.vector(), b);
  }

template <class T_FunctionSpace, class T_AdditionalFunctionSpace, class T_BilinearForm, class T_LinearForm, class T_AdditionalForm, class T_LinearSolverFactory>
void UflToNewtonLinear <T_FunctionSpace, T_AdditionalFunctionSpace, T_BilinearForm, T_LinearForm, T_AdditionalForm, T_LinearSolverFactory>::
  processAdditionalVector (dolfin::Vector& processedVec, const dolfin::Vector& vec)
  {
//TODO generalizzare al caso di MeshFunction che non sia una VertexFunction
// (se ho una linea tripla, invece di punti singoli, devo trattare interi edge, su cui vivono dof non di vertice)

/*    const dolfin::GenericDofMap& dofmap (* this->additionalFunctionSpace_.dofmap());
std::vector<dolfin::la_index> tutti_i_dof (dofmap.dofs());
for (std::size_t i=0; i!=tutti_i_dof.size(); ++i) std::cerr << tutti_i_dof[i] << ' '; std::cerr << std::endl;
    const dolfin::GenericDofMap& subdofmap (* dofmap.extract_sub_dofmap({0},*additionalFunctionSpace_.mesh()));
    std::vector<std::size_t> dofs;
    int idxs[2];
    double * vals;//[2] = {0,0};
    for (dolfin::VertexIterator v (* additionalFunctionSpace_.mesh()); !v.end(); ++v)
    {
      // if the vertex label is in labelSet_, the value is ok, hence the current iteration is skipped
      // if the vertex label is not in labelSet_, the value (0) is ok, hence the current iteration is skipped
      if (labelSet_.find((*meshFunction_)[*v]) == labelSet_.end())
        continue;
std::cerr << v->index() << ' ';

      dofmap.tabulate_entity_dofs(dofs, v->dim(), v->index());
if (dofs.size()!=2) { std::cerr << "che cavolo sta facendo??" << std::endl; exit(1);}

      // NB : GenericDofMap::tabulate_entity_dofs works with local indices, whence set_local is used, instead of set
      for (std::size_t i=0; i!=2; ++i)
      {
        idxs[i] = dofs[i];
        std::cerr << idxs[i] << ' ';
      } std::cerr << std::endl;
      vec.get_local (vals, 2, idxs);
      processedVec.set_local (vals, 2, idxs);
    }
*/
    processedVec.zero();
    for (auto it=additionalFormDofs_.begin(); it!=additionalFormDofs_.end(); ++it)
    {
      processedVec.apply("insert");
      processedVec.setitem (*it, vec[*it]);
    }
  }

} //end of namespace
#endif
