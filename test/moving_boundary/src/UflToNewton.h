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

#ifndef IVAN_UFLTONEWTON_H_INCLUDE_GUARD
#define IVAN_UFLTONEWTON_H_INCLUDE_GUARD

#include <dolfin.h>
#include <dcp/differential_problems/SubdomainType.h>
#include <dcp/differential_problems/NonlinearProblem.h>
#include "utilities.h"

//#define DAI_STAMPA
//#define ALG_STAMPA

namespace Ivan
{

  // User defined nonlinear problem
  template <class T_FunctionSpace, class T_AdditionalFunctionSpace, class T_ResidualForm, class T_JacobianForm, class T_AdditionalForm>
  class UflToNewton : public dcp::NonlinearProblem<T_ResidualForm,T_JacobianForm>,
                      public dolfin::NonlinearProblem
  {

    public:
  
      // Constructor
      UflToNewton (const dolfin::FunctionSpace& functionSpace,
                   const std::string& residualFormSolutionName,
                   const std::string& jacobianFormSolutionName = "" ) :
          dcp::NonlinearProblem<T_ResidualForm,T_JacobianForm> (functionSpace,residualFormSolutionName,jacobianFormSolutionName),
          additionalFunctionSpace_ (new T_AdditionalFunctionSpace(this->functionSpace_->mesh())),
          additionalForm_ (T_AdditionalForm (*(this->additionalFunctionSpace_))),
          bCheckFile("/u/laureandi/ifumagalli/dcp_test_output/NewtonCheck_b.pvd"),
          xCheckFile("/u/laureandi/ifumagalli/dcp_test_output/NewtonCheck_x.pvd"),
          solAlgFile("/u/laureandi/ifumagalli/dcp_test_output/solAlg.pvd"),
          solVarFile("/u/laureandi/ifumagalli/dcp_test_output/solVar.pvd"),
          timeSolAlgFile("/u/laureandi/ifumagalli/dcp_test_output/timeSolAlg.pvd"),
          timeSolVarFile("/u/laureandi/ifumagalli/dcp_test_output/timeSolVar.pvd")
      {
  /*      // Initialize class
        // Unfortunately C++ does not allow namespaces as template arguments
        init<simpleNavierStokes::FunctionSpace, simpleNavierStokes::JacobianForm,
             simpleNavierStokes::ResidualForm>(mesh, nu, gamma, dt, w);
  */
      }
  
      // Return solution function (non const version)
      // TODO sarebbe meglio tenere solo quella const di dcp::NonlinearProblem;
      //      per il momento questa serve per poter passare *solution_.vector() a dolfin::NewtonSolver
      dolfin::Function& solution ()
      {
          return this->solution_.back ().second;
      }
  
      // Overriding dcp::NonlinearProblem::setCoefficient [1]
      virtual void setCoefficient (const std::string& coefficientType, 
                                   const std::shared_ptr<const dolfin::GenericFunction> coefficientValue,
                                   const std::string& coefficientName = "default") override;

      // Overriding dcp::NonlinearProblem::setCoefficient [2]
      virtual void setCoefficient (const std::string& coefficientType,
                                   const std::shared_ptr<const dolfin::GenericFunction> coefficientValue,
                                   const std::size_t& coefficientNumber) override;

      // Overriding dcp::NonlinearProblem::setIntegrationSubdomain
      virtual void setIntegrationSubdomain (const std::string& formType,
                                             std::shared_ptr<const dolfin::MeshFunction<std::size_t>> meshFunction,
                                             const dcp::SubdomainType& subdomainType) override;
                
      // User defined residual vector (dolfin::NonlinearProblem)
      void F (dolfin::GenericVector& b, const dolfin::GenericVector& x) override;
  
      // User defined assemble of Jacobian (dolfin::NonlinearProblem)
      void J (dolfin::GenericMatrix& A, const dolfin::GenericVector& x) override;
  
      // Set solver for algebraic solution
      virtual bool setNonlinearSolver (std::shared_ptr<dolfin::NewtonSolver> solver);

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

    private:
  
  //    template<class T_FunctionSpace, class T_JacobianForm, class T_ResidualForm>
      void init (const dolfin::Mesh& mesh, const dolfin::Constant& nu,
                 const dolfin::Constant& gamma, const dolfin::Constant& dt,
                 const dolfin::Function& w)
      {
  std::cerr << "HAI USATO LA INIT!!!" << std::endl; exit(1);
  /*      // Create function space and functions
        std::shared_ptr<T_FunctionSpace> V(new T_FunctionSpace(mesh));
        solution_.reset(new Function(V));
  
        // Create forms and attach functions
        T_JacobianForm* jacobianForm_ = new T_JacobianForm(V, V);
        T_ResidualForm* residualForm_ = new T_ResidualForm(V);
        _a->u = *_u;
        _a->lmbda = lambda; _a->dt = dt; _a->theta = theta;
        _L->u = *_u; _L->u0 = *_u0;
        _L->lmbda = lambda; _L->dt = dt; _L->theta = theta;
  
        // Wrap pointers in a smart pointer
        a.reset(_a);
        L.reset(_L);
  
        // Set solution to intitial condition
        InitialConditions u_initial;
        *_u = u_initial;
  */
      }

      virtual void processAdditionalVector (dolfin::Vector& vec);
  
      // Function spaces and forms
      // (additional w.r.t. those contained in dcp::NonlinearProblem
      std::shared_ptr<dolfin::FunctionSpace> additionalFunctionSpace_;
      T_AdditionalForm additionalForm_;

      std::shared_ptr<dolfin::NewtonSolver> nonlinear_solver_;
      bool nonlinear_solver_set_;
        //TODO templatizzare/altre opzioni per avere generico solutore

dolfin::File bCheckFile;
dolfin::File xCheckFile;
dolfin::File solAlgFile;
dolfin::File solVarFile;
dolfin::File timeSolAlgFile;
dolfin::File timeSolVarFile;
  };



    // ==============================================================================================//
    // ==================================== IMPLEMENTATION ==========================================//
    // ==============================================================================================//



  template <class T_FunctionSpace, class T_AdditionalFunctionSpace, class T_ResidualForm, class T_JacobianForm, class T_AdditionalForm>
        void UflToNewton <T_FunctionSpace, T_AdditionalFunctionSpace, T_ResidualForm, T_JacobianForm, T_AdditionalForm>::
        setCoefficient (const std::string& coefficientType, 
                        const std::shared_ptr<const dolfin::GenericFunction> coefficientValue,
                        const std::string& coefficientName)
        {
            if (coefficientType == "additional_form")
            {
                dolfin::log (dolfin::DBG, "Setting jacobian form coefficient \"%s\"...", coefficientName.c_str ());
                this->additionalForm_.set_coefficient (coefficientName, coefficientValue);
            }
            else
                this->dcp::NonlinearProblem <T_ResidualForm,T_JacobianForm>::setCoefficient (coefficientType,coefficientValue,coefficientName);
        }

  template <class T_FunctionSpace, class T_AdditionalFunctionSpace, class T_ResidualForm, class T_JacobianForm, class T_AdditionalForm>
        void UflToNewton <T_FunctionSpace, T_AdditionalFunctionSpace, T_ResidualForm, T_JacobianForm, T_AdditionalForm>::
        setCoefficient (const std::string& coefficientType,
                        const std::shared_ptr<const dolfin::GenericFunction> coefficientValue,
                        const std::size_t& coefficientNumber)
        {
            if (coefficientType == "additional_form")
            {
                dolfin::log (dolfin::DBG, "Setting jacobian form coefficient number %d...", coefficientNumber);
                this->additionalForm_.set_coefficient (coefficientNumber, coefficientValue);
            }
            else
                this->dcp::NonlinearProblem <T_ResidualForm,T_JacobianForm>::setCoefficient (coefficientType,coefficientValue,coefficientNumber);
        }

  template <class T_FunctionSpace, class T_AdditionalFunctionSpace, class T_ResidualForm, class T_JacobianForm, class T_AdditionalForm>
        void UflToNewton <T_FunctionSpace, T_AdditionalFunctionSpace, T_ResidualForm, T_JacobianForm, T_AdditionalForm>::
        setIntegrationSubdomain (const std::string& formType,
                                  std::shared_ptr<const dolfin::MeshFunction<std::size_t>> meshFunction,
                                  const dcp::SubdomainType& subdomainType)
        {
            if (formType == "additional_form")
            {
                if (subdomainType == dcp::SubdomainType::INTERNAL_CELLS)
                {
                    dolfin::log (dolfin::DBG, "Setting jacobian form integration subdomain on INTERNAL_CELLS...");
                    this->additionalForm_.set_cell_domains (meshFunction);
                }
                else if (subdomainType == dcp::SubdomainType::INTERNAL_FACETS)
                {
                    dolfin::log (dolfin::DBG, "Setting jacobian form integration subdomain on INTERNAL_FACETS...");
                    this->additionalForm_.set_interior_facet_domains (meshFunction);
                }
                else if (subdomainType == dcp::SubdomainType::BOUNDARY_FACETS)
                {
                    dolfin::log (dolfin::DBG, "Setting jacobian form integration subdomain on BOUNDARY_FACETS...");
                    this->additionalForm_.set_exterior_facet_domains (meshFunction);
                }
                else if (subdomainType == dcp::SubdomainType::VERTICES)
                {
                    dolfin::log (dolfin::DBG, "Setting jacobian form integration subdomain on VERTICES...");
                    this->additionalForm_.set_vertex_domains (meshFunction);
                }
                else
                {
                    dolfin::warning ("unknown subdomain type requested while trying to apply mesh function to jacobian form"); 
                }
            }
            else
                this->dcp::NonlinearProblem <T_ResidualForm,T_JacobianForm>::setIntegrationSubdomain (formType,meshFunction,subdomainType);

        }

  template <class T_FunctionSpace, class T_AdditionalFunctionSpace, class T_ResidualForm, class T_JacobianForm, class T_AdditionalForm>
        void UflToNewton <T_FunctionSpace, T_AdditionalFunctionSpace, T_ResidualForm, T_JacobianForm, T_AdditionalForm>::
      F (dolfin::GenericVector& b, const dolfin::GenericVector& x)
      {
        // Assemble RHS (Neumann boundary conditions)
        dolfin::Assembler assembler;
        assembler.assemble (b, this->residualForm_);
std::cerr << ".h OK fino a " << __LINE__ << std::endl;

        // Take triple point into account
        dolfin::Vector vec(MPI_COMM_WORLD,b.size());
/*T_AdditionalFunctionSpace tmpFunctionSpace (this->functionSpace_->mesh());
T_AdditionalForm tmpForm (tmpFunctionSpace);
for (std::size_t i=0; i != this->additionalForm_.num_coefficients(); ++i)
{
    dolfin::Function tmpCoefficient (tmpFunctionSpace);
    tmpCoefficient = * this->additionalForm_.coefficient(i);
    tmpForm.set_coefficient (this->additionalForm_.coefficient_name(i),dolfin::reference_to_no_delete_pointer(tmpCoefficient));
}
std::cerr << ".h OK fino a " << __LINE__ << std::endl;
tmpForm.set_cell_domains (this->additionalForm_.cell_domains());
tmpForm.set_exterior_facet_domains (this->additionalForm_.exterior_facet_domains());
tmpForm.set_interior_facet_domains (this->additionalForm_.interior_facet_domains());
tmpForm.set_vertex_domains (this->additionalForm_.vertex_domains());
std::cerr << ".h OK fino a " << __LINE__ << std::endl;
//TODO sistemare: additionalFunctionSpace non subisce il mesh-moving
//     forse e' il caso di introdurre un TimeDependentFunctionSpace, cosi' tutto passa attraverso di lui...
        assembler.assemble (vec, tmpForm);*/
        assembler.assemble (vec, this->additionalForm_);
std::cerr << ".h OK fino a " << __LINE__ << std::endl;
        processAdditionalVector (vec);
#if defined(DAI_STAMPA) || defined(ALG_STAMPA)
std::cerr << "Norma residuo solo form = " << b.norm("l2") << std::endl;
#endif
        b -= vec;
//!!! controllare segno: += oppure -=

        // Apply Dirichlet boundary conditions
#if defined(DAI_STAMPA) || defined(ALG_STAMPA)
std::cerr << "Norma residuo pre-BC = " << b.norm("l2") << std::endl;
#endif
        for (auto it = this->dirichletBCs_.begin(); it != this->dirichletBCs_.end(); it++)
            {(it->second).apply (b,x); std::cerr<<"TimeBC "<<it->first<<std::endl;}
//  b.apply();
//!!! CI VUOLE (b,x) perché così metto inflow=0 nei passi interni del Newton
//!!! non faccio apply(b,x) perche' pare (a me) che quello serva se risolvo un sistema non lineare
//    tenendo come incognita la x, mentre qui, con Newton, l'incognita e' il DELTA x
//            {(it->second).apply (b, x); std::cerr<<"TimeBC "<<it->first<<std::endl;}
#if defined(DAI_STAMPA)
std::cerr << ".h OK fino a " << __LINE__ << std::endl;
dolfin::Function bFun (*this->additionalFunctionSpace_);
* bFun.vector() = b;
bCheckFile << bFun;
* bFun.vector() = x;
xCheckFile << bFun;
#endif
//dolfin::plot (bFun[0], "rhs");
/*dolfin::plot (bFun[0][0],"x rhs");
dolfin::plot (bFun[0][1],"y rhs"); dolfin::interactive();*/
#if defined(DAI_STAMPA) || defined(ALG_STAMPA)
std::cerr << " vec (size = " << vec.size() << ") = "; for (auto i=0; i!=vec.size(); ++i) std::cerr << vec[i] << ", "; std::cerr << std::endl;
std::cerr << "   b (size = " <<   b.size() << ") = "; for (auto i=0; i!=b.size(); ++i) std::cerr << b[i] << ", "; std::cerr << std::endl;
std::cerr << "Norma residuo = " << b.norm("l2") << std::endl;
std::ofstream vecFile; vecFile.open("/u/laureandi/ifumagalli/dcp_test_output/vec.csv",std::ios::out|std::ios::app);
for (auto i=0; i!=vec.size(); ++i) vecFile << vec[i] << ","; vecFile << std::endl; vecFile.close();
std::ofstream bFile; bFile.open("/u/laureandi/ifumagalli/dcp_test_output/b.csv",std::ios::out|std::ios::app);
for (auto i=0; i!=  b.size(); ++i)   bFile <<   b[i] << ",";   bFile << std::endl; bFile.close();
#endif
      }

  template <class T_FunctionSpace, class T_AdditionalFunctionSpace, class T_ResidualForm, class T_JacobianForm, class T_AdditionalForm>
        void UflToNewton <T_FunctionSpace, T_AdditionalFunctionSpace, T_ResidualForm, T_JacobianForm, T_AdditionalForm>::
      J (dolfin::GenericMatrix& A, const dolfin::GenericVector& x)
      {
        // Assemble system
        dolfin::Assembler assembler;
        assembler.assemble(A, this->jacobianForm_);
#if defined(DAI_STAMPA) || defined(ALG_STAMPA)
std::ofstream A_pre_bcFile; A_pre_bcFile.open("/u/laureandi/ifumagalli/dcp_test_output/A_pre_bc.csv",std::ios::out|std::ios::app);
for (auto i=0; i!=A.size(0); ++i) {for (auto j=0; j!=A.size(1); ++j) A_pre_bcFile << A(i,j) << ","; A_pre_bcFile << std::endl;} A_pre_bcFile << std::endl; A_pre_bcFile.close();
std::cerr << "Norma jacobiana pre-BC = " << A.norm("frobenius") << std::endl;
#endif
        for (auto it = this->dirichletBCs_.begin(); it != this->dirichletBCs_.end(); it++)
            (it->second).apply (A);
//  A.apply();
#if defined(DAI_STAMPA) || defined(ALG_STAMPA)
std::ofstream AFile; AFile.open("/u/laureandi/ifumagalli/dcp_test_output/A.csv",std::ios::out|std::ios::app);
for (auto i=0; i!=A.size(0); ++i) {for (auto j=0; j!=A.size(1); ++j) AFile << A(i,j) << ","; AFile << std::endl;} AFile << std::endl; AFile.close();
std::cerr << "   A (size = " << A.size(0) << "," << A.size(1) << ") = "; for (auto i=0; i!=A.size(0); ++i) for (auto j=0; j!=A.size(1); ++j) std::cerr << A(i,j) << ", "; std::cerr << std::endl;
std::cerr << "   A is symmetric : " << A.is_symmetric(1e-10) << std::endl;
std::cerr << "Norma jacobiana = " << A.norm("frobenius") << std::endl;
#endif
      }

  template <class T_FunctionSpace, class T_AdditionalFunctionSpace, class T_ResidualForm, class T_JacobianForm, class T_AdditionalForm>
        bool UflToNewton <T_FunctionSpace, T_AdditionalFunctionSpace, T_ResidualForm, T_JacobianForm, T_AdditionalForm>::
      setNonlinearSolver (std::shared_ptr<dolfin::NewtonSolver> solver)
        //TODO templatizzare/altre opzioni per avere generico solutore
      {
          this->nonlinear_solver_.reset(new dolfin::NewtonSolver);
          // TODO sarebbe meglio una cosa tipo la riga seguente...
          // this->nonlinear_solver_.reset (solver);
          nonlinear_solver_set_ = true;
//  nonlinear_solver_->parameters["linear_solver"] = "lu";
  nonlinear_solver_->parameters["convergence_criterion"] = "residual";
  nonlinear_solver_->parameters["maximum_iterations"] = 10;
  nonlinear_solver_->parameters["error_on_nonconvergence"] = false;
  nonlinear_solver_->parameters["relative_tolerance"] = 1e-5;
  nonlinear_solver_->parameters["absolute_tolerance"] = 1e-9;
          return nonlinear_solver_set_;
      }

  template <class T_FunctionSpace, class T_AdditionalFunctionSpace, class T_ResidualForm, class T_JacobianForm, class T_AdditionalForm>
        void UflToNewton <T_FunctionSpace, T_AdditionalFunctionSpace, T_ResidualForm, T_JacobianForm, T_AdditionalForm>::
      solve (const std::string& type)
      {
/*          std::vector<const dolfin::DirichletBC*> tmpDirichletBCs (this->dirichletBCs_.size (), nullptr);
          std::size_t counter = 0;
          for (auto i = this->dirichletBCs_.begin (); i != this->dirichletBCs_.end (); ++i)
          {
              tmpDirichletBCs[counter] = &(i->second);
              counter++;
          }

*//*          dolfin::Function solVar (* this->functionSpace_);
          * solVar.vector() = * this->solution_.back().second.vector();
      dolfin::solve (this->residualForm_==0, solVar, tmpDirichletBCs, this->jacobianForm_);//, parameters);
*//*      dolfin::solve (this->residualForm_==0, this->solution_.back().second, tmpDirichletBCs, this->jacobianForm_);//, parameters);
      return;
COSI FUNZIONA!
*/

std::cerr << "STAI USANDO L'ALGEBRICO" << std::endl;
            nonlinear_solver_->solve(*this, * this->solution_.back().second.vector());
            solAlgFile << this->solution_.back().second;
      }

  template <class T_FunctionSpace, class T_AdditionalFunctionSpace, class T_ResidualForm, class T_JacobianForm, class T_AdditionalForm>
        void UflToNewton <T_FunctionSpace, T_AdditionalFunctionSpace, T_ResidualForm, T_JacobianForm, T_AdditionalForm>::
        processAdditionalVector (dolfin::Vector& vec)
        {
            const dolfin::GenericDofMap& dofMap (* this->additionalFunctionSpace_->dofmap());
dolfin::plot(* this->additionalFunctionSpace_->mesh(), "additional mesh"); //dolfin::interactive();
dolfin::plot(* this->functionSpace_->mesh(), " non additional mesh"); //dolfin::interactive();
std::cerr << ".h OK fino a " << __LINE__ << std::endl;
            std::vector<double> coordinates (dofMap.tabulate_all_coordinates(* this->additionalFunctionSpace_->mesh()));
std::cerr << ".h OK fino a " << __LINE__ << std::endl;
//TODO sistemare: additionalFunctionSpace non subisce il mesh-moving (vedi i due plot precedenti)
//                  (a questa riga uso additionalFunctionSpace proprio perche' ha la mesh fissa, quindi posso usare
//                   SubDomains che non si muovono; il problema e' piu' su, dove c'e' un TODO simile a questo)
//     forse e' il caso di introdurre un TimeDependentFunctionSpace, cosi' tutto passa attraverso di lui...
const dolfin::GenericDofMap& dofmapComponent (*(*this->additionalFunctionSpace_)[0]->dofmap());
            std::vector<dolfin::la_index> dmCompDofs (dofmapComponent.dofs());
std::vector<dolfin::la_index> okDofs;
dolfin::Array<double> pointCoords (2);
TriplePointLeftVertex tplv;
TriplePointRightVertex tprv;
for (auto it=dmCompDofs.begin(); it!=dmCompDofs.end(); it++)
{
  pointCoords[0] = coordinates[2*(*it)];
  pointCoords[1] = coordinates[2*(*it)+1];
  if (! (tplv.inside(pointCoords,true) || tprv.inside(pointCoords,true)))
  { 
    okDofs.push_back(*it);
    const double val[1] = {0};
    const dolfin::la_index idx[1] = {*it};
    vec.set(val,1,idx);
  }
  else
    std::cerr << pointCoords[0] << "  " << pointCoords[1] << std::endl;
}
std::cerr << "neglected dofs (size " << okDofs.size() << ")" << std::endl;
std::cerr << std::endl;
        }

} //end of namespace
#endif
