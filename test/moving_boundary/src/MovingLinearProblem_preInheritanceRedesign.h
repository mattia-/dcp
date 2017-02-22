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

#ifndef IVAN_DIFFERENTIAL_PROBLEMS_MOVINGLINEARPROBLEM_H_INCLUDE_GUARD
#define IVAN_DIFFERENTIAL_PROBLEMS_MOVINGLINEARPROBLEM_H_INCLUDE_GUARD

// TODO forse in realta' questa classe potrebbe ereditare da qualcosa tipo dcp::EquationSystem

#define EVITAPLOT
#define ADDITIONALVECTOR

#include <dolfin.h>
#include <dcp/differential_problems/SubdomainType.h>
#include <dcp/differential_problems/LinearProblem.h>
#include <fstream>
//#include "utilities.h"
#include "MeshManager.h"

extern struct ProblemData problemData;

namespace Ivan
{

    /*! \class MovingLinearProblem MovingLinearProblem.h
     *  \brief Class for linear differential problems with moving domain.
     *
     *  This class represents problem of the form
     *  \f[
     *      \mbox{Find } u \in V \left(\Omega\right) : 
     *      a_\Omega \left(u , v\right) 
     *      = 
     *      F_\Omega \left(v\right) 
	 *		+
	 *		\widetilde{F}_{\widetilde{\Omega}} \left(v\right)
     *      \ \forall\,v\,\in\,V \left(\Omega\right)
     *  \f]
     *  with \f$ \Omega, \widtilde{\Omega} \f$ two domains in \f$ \mathbb{R}^d \f$,
	 *	\f$ V \left(\Omega\right)\f$ Hilbert spaces of functions over the domain \f$\Omega\f$,
	 *  \f$ a_\Omega \left(u , v\right) : V \left(\Omega\right) \times V \left(\Omega\right) \rightarrow \mathds{R}\f$ generic bilinear form on \f$V \left(\Omega\right)\f$, with integration occurring over \f$ \Omega \f$,
     *  and \f$ F_\Omega \left(v\right) : V \rightarrow \mathds{R}, \widetilde{F}_{\widetilde{\Omega} \left(v\right) : V \rightarrow \mathds{R} \f$ linear forms on the same space, with integration occurring over \f$ \Omega , \widetilde{\Omega} \f$, respectively.
     *  
     *  It inherits publicly from \c LinearProblem and it extends
     *  its functionalities to a moving-domain differential problem.
     */

template <class T_FunctionSpace, class T_AdditionalFunctionSpace, class T_BilinearForm, class T_LinearForm, class T_PreviousForm, class T_AdditionalForm, class T_LinearSolverFactory = dcp::LinearSolverFactory>
class MovingLinearProblem : public dcp::LinearProblem <T_BilinearForm, T_LinearForm, T_LinearSolverFactory>
{

  public:

    typedef T_FunctionSpace FunctionSpace;
    typedef T_AdditionalFunctionSpace AdditionalFunctionSpace;
    typedef T_BilinearForm BilinearForm;
    typedef T_LinearForm LinearForm;
    typedef T_PreviousForm PreviousForm;
    typedef T_AdditionalForm AdditionalForm;
    typedef T_LinearSolverFactory LinearSolverFactory;

    // Constructor
    MovingLinearProblem (std::shared_ptr<dolfin::FunctionSpace> functionSpace, const std::vector<dolfin::la_index> & additionalFormDofs);

    // Return solution function (non const version)
    // TODO sarebbe meglio tenere solo quella const di dcp::LinearProblem;
    //      per il momento questa serve per poter passare *solution_.vector() a dolfin::LinearSolver
    dolfin::Function& solution ()
		{ return this->solution_.back ().second; }


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
	// TODO: al momento e' praticamente inutile... vedi se c'e' una policy simile in altri figli di AbstractProblem
    virtual bool setLinearSolver (std::shared_ptr<dolfin::GenericLinearSolver> solver);

	//! Assemble parts of the problem depending on the previous mesh
	/*!
	 *	This method assembles \c previousForm_ and stores it in \c preassembled_
	 */
    virtual void preassemble ();

    //! Solve problem
    /*!
     *  This method solves the problem defined. It uses the private members' value to set the problem and then
     *  stores the solution in the private member \c solution_. See documentation of \c dcp::AbstractProblem
     *  for more details on how the protected member \c solution_ works and why it is declared as a 
     *  \c std::pair.
     *
     *  \param type the solution type requested. In this class, the only possibilities are to set
     *  \c type equal to \c "default" or "save_residual". In the latter case, the
  	 * 	the residual vector of the system is stored in residual.dat
     */
    virtual void solve (const std::string& type = "default") override;

    //! Set mesh-and-dofs manager
    virtual void setMeshManager (const MeshManager<dolfin::ALE> & meshManager)
		{ meshManager_ = & meshManager; }

    //! Update the problem after mesh changing
    /*! After the mesh has changed, the forms, function spaces and functions defined on it have to be updated.
     *  This update is performed via dolfin::adapt
     */
    virtual void adapt ();
    //TODO remove! if not necessary

    T_AdditionalForm & additionalForm ()
		{ return additionalForm_; }

    const std::vector<dolfin::la_index> & additionalFormDofs ()
    { return additionalFormDofs_; }

  protected:

	//! Store in processedVec only the interesting dofs of \c inputVec
	/*!
	 *	This method sets \c processedVec as zero everywhere but in the entries
	 * 	identified by \c additionalFormDofs_, where it is set
	 *	equal to \c inputVec.
	 */
    virtual void processAdditionalVector (dolfin::Vector& processedVec, const dolfin::Vector& inputVec);

    const std::vector<dolfin::la_index> & additionalFormDofs_;
    dolfin::FunctionSpace additionalFunctionSpace_;
    T_AdditionalForm additionalForm_;

    std::shared_ptr<dolfin::GenericLinearSolver> linear_solver_;
    bool linear_solver_set_;

    T_PreviousForm previousForm_;
    dolfin::Vector preassembled_;
  
    const MeshManager<dolfin::ALE> * meshManager_;

}; //end of class definition


    // ==============================================================================================//
    // ==================================== IMPLEMENTATION ==========================================//
    // ==============================================================================================//

template <class T_FunctionSpace, class T_AdditionalFunctionSpace, class T_BilinearForm, class T_LinearForm, class T_PreviousForm, class T_AdditionalForm, class T_LinearSolverFactory>
MovingLinearProblem <T_FunctionSpace, T_AdditionalFunctionSpace, T_BilinearForm, T_LinearForm, T_PreviousForm, T_AdditionalForm, T_LinearSolverFactory>::
  MovingLinearProblem (std::shared_ptr<dolfin::FunctionSpace> functionSpace, const std::vector<dolfin::la_index> & additionalFormDofs) :
    dcp::LinearProblem<T_BilinearForm,T_LinearForm,T_LinearSolverFactory> (functionSpace),
    additionalFormDofs_ (additionalFormDofs),
    additionalFunctionSpace_ (T_AdditionalFunctionSpace(this->functionSpace_->mesh())),
		// In this way, the mesh object is shared between \c additionalFunctionSpace_
		// and \c functionSpace_, so that moving the mesh affects both.
    additionalForm_ (T_AdditionalForm (this->additionalFunctionSpace_)),
    previousForm_ (T_PreviousForm (functionSpace))
  {
          std::ofstream outTP (problemData.savepath+"tpValues.csv", std::ofstream::out);
          outTP << "tpVal, \n0, " << std::endl;
          outTP.flush();
          outTP.close();
          std::ifstream outTPread (problemData.savepath+"tpValues.csv", std::ifstream::in);
          outTPread.get();
          std::ofstream residualFile (problemData.savepath+"residual.dat",std::ofstream::out);
          residualFile.close();
  }

template <class T_FunctionSpace, class T_AdditionalFunctionSpace, class T_BilinearForm, class T_LinearForm, class T_PreviousForm, class T_AdditionalForm, class T_LinearSolverFactory>
void MovingLinearProblem <T_FunctionSpace, T_AdditionalFunctionSpace, T_BilinearForm, T_LinearForm, T_PreviousForm, T_AdditionalForm, T_LinearSolverFactory>::
  setCoefficient (const std::string& coefficientType, 
                  const std::shared_ptr<const dolfin::GenericFunction> coefficientValue,
                  const std::string& coefficientName)
  {
      if (coefficientType == "additional_form")
      {
          dolfin::log (dolfin::DBG, "Setting linear form coefficient \"%s\"...", coefficientName.c_str ());
          this->additionalForm_.set_coefficient (coefficientName, coefficientValue);
      }
      else if (coefficientType == "previous_form")
      {
          dolfin::log (dolfin::DBG, "Setting linear form coefficient \"%s\"...", coefficientName.c_str ());
          this->previousForm_.set_coefficient (coefficientName, coefficientValue);
      }
      else
          this->dcp::LinearProblem <T_BilinearForm,T_LinearForm,T_LinearSolverFactory>::setCoefficient (coefficientType,coefficientValue,coefficientName);
  }

template <class T_FunctionSpace, class T_AdditionalFunctionSpace, class T_BilinearForm, class T_LinearForm, class T_PreviousForm, class T_AdditionalForm, class T_LinearSolverFactory>
void MovingLinearProblem <T_FunctionSpace, T_AdditionalFunctionSpace, T_BilinearForm, T_LinearForm, T_PreviousForm, T_AdditionalForm, T_LinearSolverFactory>::
  setCoefficient (const std::string& coefficientType,
                  const std::shared_ptr<const dolfin::GenericFunction> coefficientValue,
                  const std::size_t& coefficientNumber)
  {
      if (coefficientType == "additional_form")
      {
          dolfin::log (dolfin::DBG, "Setting linear form coefficient number %d...", coefficientNumber);
          this->additionalForm_.set_coefficient (coefficientNumber, coefficientValue);
      }
      else if (coefficientType == "previous_form")
      {
          dolfin::log (dolfin::DBG, "Setting linear form coefficient number %d...", coefficientNumber);
          this->previousForm_.set_coefficient (coefficientNumber, coefficientValue);
      }
      else
          this->dcp::LinearProblem <T_BilinearForm,T_LinearForm,T_LinearSolverFactory>::setCoefficient (coefficientType,coefficientValue,coefficientNumber);
  }

template <class T_FunctionSpace, class T_AdditionalFunctionSpace, class T_BilinearForm, class T_LinearForm, class T_PreviousForm, class T_AdditionalForm, class T_LinearSolverFactory>
void MovingLinearProblem <T_FunctionSpace, T_AdditionalFunctionSpace, T_BilinearForm, T_LinearForm, T_PreviousForm, T_AdditionalForm, T_LinearSolverFactory>::
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
      else if (formType == "previous_form")
      {
          if (subdomainType == dcp::SubdomainType::INTERNAL_CELLS)
          {
              dolfin::log (dolfin::DBG, "Setting linear form integration subdomain on INTERNAL_CELLS...");
              this->previousForm_.set_cell_domains (meshFunction);
          }
          else if (subdomainType == dcp::SubdomainType::INTERNAL_FACETS)
          {
              dolfin::log (dolfin::DBG, "Setting linear form integration subdomain on INTERNAL_FACETS...");
              this->previousForm_.set_interior_facet_domains (meshFunction);
          }
          else if (subdomainType == dcp::SubdomainType::BOUNDARY_FACETS)
          {
              dolfin::log (dolfin::DBG, "Setting linear form integration subdomain on BOUNDARY_FACETS...");
              this->previousForm_.set_exterior_facet_domains (meshFunction);
          }
          else if (subdomainType == dcp::SubdomainType::VERTICES)
          {
              dolfin::log (dolfin::DBG, "Setting linear form integration subdomain on VERTICES...");
              this->previousForm_.set_vertex_domains (meshFunction);
          }
          else
          {
              dolfin::warning ("unknown subdomain type requested while trying to apply mesh function to linear form"); 
          }
      }
      else
          this->dcp::LinearProblem <T_BilinearForm,T_LinearForm,T_LinearSolverFactory>::setIntegrationSubdomain (formType,meshFunction,subdomainType);
  }

template <class T_FunctionSpace, class T_AdditionalFunctionSpace, class T_BilinearForm, class T_LinearForm, class T_PreviousForm, class T_AdditionalForm, class T_LinearSolverFactory>
bool MovingLinearProblem <T_FunctionSpace, T_AdditionalFunctionSpace, T_BilinearForm, T_LinearForm, T_PreviousForm, T_AdditionalForm, T_LinearSolverFactory>::
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

template <class T_FunctionSpace, class T_AdditionalFunctionSpace, class T_BilinearForm, class T_LinearForm, class T_PreviousForm, class T_AdditionalForm, class T_LinearSolverFactory>
void MovingLinearProblem <T_FunctionSpace, T_AdditionalFunctionSpace, T_BilinearForm, T_LinearForm, T_PreviousForm, T_AdditionalForm, T_LinearSolverFactory>::
  preassemble ()
  {
const dolfin::GenericVector & oldVec ((this->solution_.size()>1) ? (* this->solution_[this->solution_.size()-2].second.vector()) : (* this->solution_.back().second.vector()));
/*std::cerr << "preassVec (" << oldVec.size() << ") : ";
for (std::size_t i (0); i!=oldVec.size(); ++i)
  std::cerr << oldVec[i] << '=' << oldVec[i] << ", ";
std::cerr << std::endl;*/
        dolfin::Assembler assembler;

        assembler.assemble (this->preassembled_, this->previousForm_);
  }

template <class T_FunctionSpace, class T_AdditionalFunctionSpace, class T_BilinearForm, class T_LinearForm, class T_PreviousForm, class T_AdditionalForm, class T_LinearSolverFactory>
void MovingLinearProblem <T_FunctionSpace, T_AdditionalFunctionSpace, T_BilinearForm, T_LinearForm, T_PreviousForm, T_AdditionalForm, T_LinearSolverFactory>::
  solve (const std::string& type)
  {
        dolfin::Assembler assembler;

        dolfin::Matrix A;
        assembler.assemble (A, this->bilinearForm_);

        dolfin::Vector b;
        assembler.assemble (b, this->linearForm_);

        b += this->preassembled_;

        dolfin::Vector vec (MPI_COMM_WORLD,b.size());
#ifndef EVITAPLOT
  dolfin::plot (*additionalFunctionSpace_.mesh(),"additional mesh"); //dolfin::interactive ();
#endif
        assembler.assemble (vec, this->additionalForm_);
        dolfin::Vector processedVec (MPI_COMM_WORLD,vec.size());
        processedVec.zero();
#ifdef ADDITIONALVECTOR
        processAdditionalVector (processedVec, vec);
#endif
        b += processedVec;

        for (auto it = this->dirichletBCs_.begin(); it != this->dirichletBCs_.end(); it++)
        {
          (it->second).apply (A,b);
//          std::cerr<<"TimeBC "<<it->first<<std::endl;
        }

        linear_solver_->solve(A, * this->solution_.back().second.vector(), b);

        if (type == "save_residual")
        {
          dolfin::Vector residual (MPI_COMM_WORLD,b.size());
          A.mult (* this->solution_.back().second.vector(), residual);
          residual -= b;
          double resdata[residual.size()];// (residual.data());
          std::vector<dolfin::la_index> idxs (residual.size());
          std::iota (idxs.begin(), idxs.end(), 0);
          residual.get (resdata, residual.size(), &idxs[0]);
          std::ofstream residualFile (problemData.savepath+"residual.dat",std::ofstream::app);
          for (std::size_t i (0); i<residual.size(); ++i)
            residualFile << resdata[i] << ',';
          residualFile << std::endl;
          residualFile.close();
        }

        std::ofstream outTP (problemData.savepath+"tpValues.csv", std::ofstream::app);
        outTP << processedVec.inner (* this->solution_.back().second.vector()) / problemData.dt << ", " << std::endl;
        outTP.close();
  }

template <class T_FunctionSpace, class T_AdditionalFunctionSpace, class T_BilinearForm, class T_LinearForm, class T_PreviousForm, class T_AdditionalForm, class T_LinearSolverFactory>
void MovingLinearProblem <T_FunctionSpace, T_AdditionalFunctionSpace, T_BilinearForm, T_LinearForm, T_PreviousForm, T_AdditionalForm, T_LinearSolverFactory>::
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
/*    for (auto it=additionalFormDofs_.begin(); it!=additionalFormDofs_.end(); ++it)
    {
      processedVec.apply("insert");
      processedVec.setitem (*it, vec[*it]);
std::cerr << "processed values " << vec[*it] << ", ";
    }
*/
    processedVec.setitem (additionalFormDofs_.back(), (90==problemData.thetaS ? 0 : problemData.dt*problemData.gamma*cos(problemData.thetaS*3.14159265/180.0)*problemData.lx));
      // TODO generalizzazione geometrica (secondo l'idea del for di sopra,
      //      che pero' da' dei nan se si raffina troppo in spazio (forse per il /h
      //      che c'e' in computeFreeSurfaceStress_onlyTP.ufl))
      // Al momento questa riga vale perche' il muro e' verticale (=> v.t = v_y) e perche' siamo in cilindriche (*lx)
std::cerr << "processed values " << processedVec[additionalFormDofs_.back()] << ", ";
std::cerr << std::endl;
  }

template <class T_FunctionSpace, class T_AdditionalFunctionSpace, class T_BilinearForm, class T_LinearForm, class T_PreviousForm, class T_AdditionalForm, class T_LinearSolverFactory>
void MovingLinearProblem <T_FunctionSpace, T_AdditionalFunctionSpace, T_BilinearForm, T_LinearForm, T_PreviousForm, T_AdditionalForm, T_LinearSolverFactory>::
  adapt ()
  {
/*    this->bilinearForm_ = static_cast<const T_BilinearForm &> (dolfin::adapt (this->bilinearForm_, meshManager_->mesh()));
    this->linearForm_ = static_cast<const T_LinearForm &> (dolfin::adapt (this->linearForm_, meshManager_->mesh()));
    this->additionalForm_ = static_cast<const T_AdditionalForm &> (dolfin::adapt (this->additionalForm_, meshManager_->mesh()));
    this->previousForm_ = static_cast<const T_PreviousForm &> (dolfin::adapt (this->previousForm_, meshManager_->mesh()));*/
std::cerr << " OK fino a " << __LINE__ << std::endl;
    const std::shared_ptr<const dolfin::Mesh> mesh (dolfin::reference_to_no_delete_pointer (* this->functionSpace_->mesh()));
std::cerr << " OK fino a " << __LINE__ << std::endl;
    dolfin::adapt (this->bilinearForm_, mesh);
std::cerr << " OK fino a " << __LINE__ << std::endl;
    dolfin::adapt (this->linearForm_, mesh);
std::cerr << " OK fino a " << __LINE__ << std::endl;
    dolfin::adapt (this->additionalForm_, mesh);
std::cerr << " OK fino a " << __LINE__ << std::endl;
    dolfin::adapt (this->previousForm_, mesh);
std::cerr << " OK fino a " << __LINE__ << std::endl;
std::cerr << "Si', ha usato Ivan::MovingLinearProblem<...>::adapt()" << std::endl;
    exit(10);
  }

} //end of namespace
#endif
