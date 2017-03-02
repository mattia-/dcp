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

#ifndef DCP_DIFFERENTIAL_PROBLEMS_MOVINGLINEARPROBLEM_H_INCLUDE_GUARD
#define DCP_DIFFERENTIAL_PROBLEMS_MOVINGLINEARPROBLEM_H_INCLUDE_GUARD

// TODO forse in realta' questa classe potrebbe ereditare da qualcosa tipo dcp::EquationSystem

#define EVITAPLOT
#define ADDITIONALVECTOR

#include <dolfin.h>
#include <dcp/differential_problems/SubdomainType.h>
#include <dcp/differential_problems/LinearProblem.h>
#include <fstream>
#include <dcp/differential_problems/MovingAbstractProblem.h>
#include <dcp/differential_problems/AdditionalProcessor.h>

extern struct ProblemData problemData;

namespace dcp
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
     *  It inherits publicly from #LinearProblem AND #MovingAbstractProblem in order to
     *  extend the functionalities of #LinearProblem to the moving-domain case.
     */

template <class T_FunctionSpace, class T_AdditionalFunctionSpace, class T_BilinearForm, class T_LinearForm, class T_PreviousForm, class T_AdditionalForm, class T_LinearSolverFactory = dcp::LinearSolverFactory>
class MovingLinearProblem : public dcp::LinearProblem <T_BilinearForm, T_LinearForm, T_LinearSolverFactory>,
                            public dcp::MovingAbstractProblem
{

  public:

    typedef T_FunctionSpace FunctionSpace;
    typedef T_AdditionalFunctionSpace AdditionalFunctionSpace;
    typedef T_BilinearForm BilinearForm;
    typedef T_LinearForm LinearForm;
    typedef T_PreviousForm PreviousForm;
    typedef T_AdditionalForm AdditionalForm;
    typedef T_LinearSolverFactory LinearSolverFactory;

    //! Constructor
    /*!
     * @param functionSpace The problem FE space (@see #LinearProblem)
     * @param additionalFormDofs Dofs for additional terms, to be added to the rhs of the problem
     * @param additionalProcessor (passed through std::shared_ptr) Functor processing the values coming from the assembling of a \c T_AdditionalForm, before adding them to the rhs of the problem (@see #AdditionalProcessor)
     */
    MovingLinearProblem (std::shared_ptr<dolfin::FunctionSpace> functionSpace, const std::vector<dolfin::la_index> & additionalFormDofs, std::shared_ptr<AdditionalProcessor> additionalProcessor);

    //! Auxiliary constructor for the case when no additional terms must be added to the algebraic problem
    /*!
     * It behaves exactly like passing an empty vector as the additionalFormDofs parameter and a do-nothing additionalProcessor in the other ctor.
     * TODO: this issue should be taken into account also at the template-arguments level
     */
    MovingLinearProblem (std::shared_ptr<dolfin::FunctionSpace> functionSpace);

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
    virtual void preassemble () override;

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
    virtual void setMeshManager (const dcp::MeshManager<dolfin::ALE> & meshManager)
		{ meshManager_ = & meshManager; }

    //! Update the problem after mesh changing
    /*! After the mesh has changed, the forms, function spaces and functions defined on it have to be updated.
     *  This update is performed via dolfin::adapt
     */
    virtual void adapt ();
    //TODO remove! if not necessary

    //T_AdditionalForm & additionalForm ()
    dolfin::Form & additionalForm ()
		{ return additionalForm_; }

    const std::vector<dolfin::la_index> & additionalFormDofs ()
    { return additionalFormDofs_; }

    //! TODO Clone method. Overrides method in \c AbstractProblem
    /*!
     *  Still to be implemented.
     *
     *  It uses the parameter \c clone_method to decide which type of cloning to perform.
     *  Possible values for such parameter are:
     *  \li deep_clone the new object is created calling the constructor that takes a mesh and a function 
     *  space as input, thus creating a copy of such objects and returning a completely independent 
     *  cloned object. 
     *  \li shallow_clone calls the constructor that takes shared pointers as input: the mesh and
     *  the function space are not copied but shared between the current object and its clone. 
     *  
     *  The default value for parameter \c clone_method is \c shallow_clone
     *  
     *  \return a pointer to the cloned object
     */
    virtual dcp::MovingLinearProblem<T_FunctionSpace,T_AdditionalFunctionSpace,T_BilinearForm,T_LinearForm,T_PreviousForm,T_AdditionalForm,T_LinearSolverFactory>* clone () const override;

  protected:

	//! Store in processedVec only the interesting dofs of \c inputVec
	/*!
	 *	This method sets \c processedVec as zero everywhere but in the entries
	 * 	identified by \c additionalFormDofs_, where it is set
	 *	by means of #additionalProcessor_.
	 */
    virtual void processAdditionalVector (dolfin::Vector& processedVec, const dolfin::Vector& inputVec);

    const std::vector<dolfin::la_index> & additionalFormDofs_;
//AFS//    dolfin::FunctionSpace additionalFunctionSpace_;
    //T_AdditionalForm additionalForm_;
    std::shared_ptr<T_AdditionalForm> additionalForm_PUNTATORECHEMISERVEPERFAREDELLEPROVEPRIMADIPASSAREadditionalForm_DAFUORINELCOSTRUTTORE;
    dolfin::Form & additionalForm_;
    std::shared_ptr<AdditionalProcessor> additionalProcessor_;

    std::shared_ptr<dolfin::GenericLinearSolver> linear_solver_;
    bool linear_solver_set_;

    T_PreviousForm previousForm_;
    dolfin::Vector preassembled_;
  
    const dcp::MeshManager<dolfin::ALE> * meshManager_;

}; //end of class definition


    // ==============================================================================================//
    // ==================================== IMPLEMENTATION ==========================================//
    // ==============================================================================================//

template <class T_FunctionSpace, class T_AdditionalFunctionSpace, class T_BilinearForm, class T_LinearForm, class T_PreviousForm, class T_AdditionalForm, class T_LinearSolverFactory>
MovingLinearProblem <T_FunctionSpace, T_AdditionalFunctionSpace, T_BilinearForm, T_LinearForm, T_PreviousForm, T_AdditionalForm, T_LinearSolverFactory>::
  MovingLinearProblem (std::shared_ptr<dolfin::FunctionSpace> functionSpace) :
    MovingLinearProblem<T_FunctionSpace,T_AdditionalFunctionSpace,T_BilinearForm,T_LinearForm,T_PreviousForm,T_AdditionalForm,T_LinearSolverFactory> (functionSpace, std::vector<dolfin::la_index>(), std::make_shared<AdditionalProcessor>(AdditionalProcessor()))
  {
  }
template <class T_FunctionSpace, class T_AdditionalFunctionSpace, class T_BilinearForm, class T_LinearForm, class T_PreviousForm, class T_AdditionalForm, class T_LinearSolverFactory>
MovingLinearProblem <T_FunctionSpace, T_AdditionalFunctionSpace, T_BilinearForm, T_LinearForm, T_PreviousForm, T_AdditionalForm, T_LinearSolverFactory>::
  MovingLinearProblem (std::shared_ptr<dolfin::FunctionSpace> functionSpace, const std::vector<dolfin::la_index> & additionalFormDofs, std::shared_ptr<AdditionalProcessor> additionalProcessor) :
    dcp::AbstractProblem(functionSpace),
      // explicit call to "grandparent" ctor req'd by virtual inheritance
    dcp::MovingAbstractProblem(functionSpace),
      // explicit call to ALL parents ctor req'd by virtual inheritance
    dcp::LinearProblem<T_BilinearForm,T_LinearForm,T_LinearSolverFactory> (functionSpace),
    additionalFormDofs_ (additionalFormDofs),
//AFS//    additionalFunctionSpace_ (T_AdditionalFunctionSpace(this->functionSpace_->mesh())),
		// In this way, the mesh object is shared between \c additionalFunctionSpace_
		// and \c functionSpace_, so that moving the mesh affects both.
//AFS//    additionalForm_PUNTATORECHEMISERVEPERFAREDELLEPROVEPRIMADIPASSAREadditionalForm_DAFUORINELCOSTRUTTORE (new T_AdditionalForm (this->additionalFunctionSpace_)),
    additionalForm_PUNTATORECHEMISERVEPERFAREDELLEPROVEPRIMADIPASSAREadditionalForm_DAFUORINELCOSTRUTTORE (new T_AdditionalForm (std::shared_ptr<dolfin::FunctionSpace> (new T_AdditionalFunctionSpace (this->functionSpace_->mesh())))),
    additionalForm_ (* additionalForm_PUNTATORECHEMISERVEPERFAREDELLEPROVEPRIMADIPASSAREadditionalForm_DAFUORINELCOSTRUTTORE), //T_AdditionalForm (this->additionalFunctionSpace_)),
    additionalProcessor_ (additionalProcessor),
    linear_solver_ (new dolfin::LinearSolver),
    linear_solver_set_ (false),
    previousForm_ (T_PreviousForm (functionSpace))
  {
          std::ofstream outAddVal (problemData.savepath+"additionalValues.csv", std::ofstream::out);
          outAddVal << "additionalValues, \n0, " << std::endl;
          outAddVal.flush();
          outAddVal.close();
          std::ifstream outAddValread (problemData.savepath+"additionalValues.csv", std::ifstream::in);
          outAddValread.get();
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
//AFS//  dolfin::plot (*additionalFunctionSpace_.mesh(),"additional mesh"); //dolfin::interactive ();
#endif
        assembler.assemble (vec, this->additionalForm_);
        dolfin::Vector processedVec (MPI_COMM_WORLD,vec.size());
        processedVec.zero();
#ifdef ADDITIONALVECTOR
        processAdditionalVector (processedVec, vec);
#endif
        b += processedVec;

        for (auto it = this->dirichletBCs_.begin(); it != this->dirichletBCs_.end(); it++)
          (it->second).apply (A,b);

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

        std::ofstream outAddVal (problemData.savepath+"additionalValues.csv", std::ofstream::app);
        outAddVal << processedVec.inner (* this->solution_.back().second.vector()) / problemData.dt << ", " << std::endl;
        outAddVal.close();
  }

template <class T_FunctionSpace, class T_AdditionalFunctionSpace, class T_BilinearForm, class T_LinearForm, class T_PreviousForm, class T_AdditionalForm, class T_LinearSolverFactory>
void MovingLinearProblem <T_FunctionSpace, T_AdditionalFunctionSpace, T_BilinearForm, T_LinearForm, T_PreviousForm, T_AdditionalForm, T_LinearSolverFactory>::
  processAdditionalVector (dolfin::Vector& processedVec, const dolfin::Vector& vec)
  {
//TODO generalizzare al caso di MeshFunction che non sia una VertexFunction
// (se ho una linea tripla, invece di punti singoli, devo trattare interi edge, su cui vivono dof non di vertice)

    processedVec.zero();

    (* additionalProcessor_) (processedVec, vec, additionalFormDofs_);
  }

template <class T_FunctionSpace, class T_AdditionalFunctionSpace, class T_BilinearForm, class T_LinearForm, class T_PreviousForm, class T_AdditionalForm, class T_LinearSolverFactory>
void MovingLinearProblem <T_FunctionSpace, T_AdditionalFunctionSpace, T_BilinearForm, T_LinearForm, T_PreviousForm, T_AdditionalForm, T_LinearSolverFactory>::
  adapt ()
  {
    const std::shared_ptr<const dolfin::Mesh> mesh (dolfin::reference_to_no_delete_pointer (* this->functionSpace_->mesh()));
    dolfin::adapt (this->bilinearForm_, mesh);
    dolfin::adapt (this->linearForm_, mesh);
    dolfin::adapt (this->additionalForm_, mesh);
    dolfin::adapt (this->previousForm_, mesh);
std::cerr << "Si', ha usato dcp::MovingLinearProblem<...>::adapt()" << std::endl;
    exit(10);
  }

template <class T_FunctionSpace, class T_AdditionalFunctionSpace, class T_BilinearForm, class T_LinearForm, class T_PreviousForm, class T_AdditionalForm, class T_LinearSolverFactory>
MovingLinearProblem<T_FunctionSpace,T_AdditionalFunctionSpace,T_BilinearForm,T_LinearForm,T_PreviousForm,T_AdditionalForm,T_LinearSolverFactory>* MovingLinearProblem <T_FunctionSpace, T_AdditionalFunctionSpace, T_BilinearForm, T_LinearForm, T_PreviousForm, T_AdditionalForm, T_LinearSolverFactory>::
  clone () const
  {
    std::cerr << std::endl << "dcp::MovingLinearProblem::clone() is TODO" << std::endl;
    exit(10);
  }

} //end of namespace
#endif
