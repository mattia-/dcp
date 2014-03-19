#ifndef SRC_DIFFERENTIALPROBLEM_NONLINEARDIFFERENTIALPROBLEM_HPP_INCLUDE_GUARD
#define SRC_DIFFERENTIALPROBLEM_NONLINEARDIFFERENTIALPROBLEM_HPP_INCLUDE_GUARD

#include <dolfin.h>
#include <vector>
#include <map>
#include <string>
#include <memory>
#include <DifferentialProblem/AbstractDifferentialProblem.hpp>

namespace control_problem
{
	/*! \class NonLinearDifferentialProblem NonLinearDifferentialProblem.hpp
	 *  \brief Class for non-linear differential problems.
	 *
	 *  This class represents problem of the form
	 *  \f[
	 *      \mbox{Find } u \in V : F \left(u, v\right) = 0 \ \forall\,v\,\in\,V
	 *  \f]
	 *  with \f$ F \left(u, v\right) : V \times V \rightarrow \mathds{R}\f$ non linear in the problem unknown \f$u\f$.
	 *  
	 *  It inherits publicly from \c AbstractDifferentialProblem
	 *  and it extends its functionalities to a concrete differential
	 *  problem.
	 *  Template arguments are:
	 *  \arg T_ResidualForm the residual form type, that is the type that describes the form \f$F\f$
	 *  \arg T_JacobianForm the jacobian form type, that is the type of the derivative of \f$F\f$ with respect to \f$u\f$,
	 *  used in the Newton-Raphson method iterations
	 *  This is needed because the dolfin function "coefficient_number" was neither declared
	 *  virtual nor implemented by FEniCS's developers in the dolfin class "Form", from which
	 *  we first thought to derive for the protected members of the NonLinearDifferentialProblem
	 *  class
	 */

	template <class T_ResidualForm, class T_JacobianForm>
		class NonLinearDifferentialProblem : public AbstractDifferentialProblem
		{
			// ---------------------------------------------------------------------------------------------//	

			public:
				/******************* CONSTRUCTORS *******************/
				//! Default constructor is deleted. The class is not default constructable.
				NonLinearDifferentialProblem () = delete;

				//!  Constructor with shared pointers
				/*!
			 	 *  \param mesh the problem mesh as a const std::shared_ptr to dolfin::Mesh
			 	 *  \param functionSpace the problem finite element space as a const std::shared_ptr to dolfin::FunctionSpace
				 *  \param residualFormSolutionName a string that identifies the name of the function representing the problem
				 *  solution in the ResidualForm 
				 *  \param jacobianFormSolutionName a string that identifies the name of the function representing the problem
				 *  solution in the JacobianForm. Default value is the empty string, in which case \c residualFormSolutionName_
				 *  will be used for jacobianFormSolutionName_
			 	 *  The stored mesh's and function space's ownership will be shared between the object and the input argument.
			 	 *  The bilinear and linear form will be created too, calling the constructor which takes the function space
			 	 *  as input.
			 	 */
				NonLinearDifferentialProblem (const std::shared_ptr<dolfin::Mesh> mesh, 
				                              const std::shared_ptr<dolfin::FunctionSpace> functionSpace,
			                           	      const std::string& residualFormSolutionName,
			                           	      const std::string& jacobianFormSolutionName);
				

				//! Constructor with references
				/*!
			 	 *  \param mesh the problem mesh as a const dolfin::Mesh&
			 	 *  \param functionSpace the problem finite element space as a const dolfin::FunctionSpace&
				 *  \param residualFormSolutionName a string that identifies the name of the function representing the problem
				 *  solution in the ResidualForm 
				 *  \param jacobianFormSolutionName a string that identifies the name of the function representing the problem
				 *  solution in the JacobianForm. Default value is the empty string, in which case \c residualFormSolutionName_
				 *  will be used for jacobianFormSolutionName_
			 	 *  The stored mesh's and function space's ownership will be unique to the object, since the pointers are 
			 	 *  initialized using the \c new operator and mesh's and functionSpace's copy constructor
			 	 *  The bilinear and linear form will be created too, calling the constructor which takes the function space
			 	 *  as input.
			 	 */
				NonLinearDifferentialProblem (const dolfin::Mesh& mesh, 
			                           	      const dolfin::FunctionSpace& functionSpace,
			                           	      const std::string& residualFormSolutionName,
			                           	      const std::string& jacobianFormSolutionName);

				//! Constructor with rvalue references
				/*!
			 	 *  \param mesh the problem mesh as a dolfin::Mesh&&
			 	 *  \param functionSpace the problem finite element space as a dolfin::FunctionSpace&&
				 *  \param residualFormSolutionName a string that identifies the name of the function representing the problem
				 *  solution in the ResidualForm 
				 *  \param jacobianFormSolutionName a string that identifies the name of the function representing the problem
				 *  solution in the JacobianForm. Default value is the empty string, in which case \c residualFormSolutionName_
				 *  will be used for jacobianFormSolutionName_
			 	 *  The stored mesh's and function space's ownership will be unique to the object, since the pointers are 
			 	 *  initialized using the \c new operator and mesh's and functionSpace's move constructor
			 	 *  The bilinear and linear form will be created too, calling the constructor which takes the function space
			 	 *  as input.
			 	 */
				NonLinearDifferentialProblem (dolfin::Mesh&& mesh, 
			                           	      dolfin::FunctionSpace&& functionSpace,
			                           	      const std::string& residualFormSolutionName,
			                           	      const std::string& jacobianFormSolutionName);

				
				//!  Constructor with shared pointers [2]
				/*!
			 	 *  \param mesh the problem mesh as a const std::shared_ptr to dolfin::Mesh
			 	 *  \param functionSpace the problem finite element space as a const std::shared_ptr to dolfin::FunctionSpace
				 *  \param residualForm a const reference to the problem's bilinear form
				 *  \param jacobianForm a const reference to the problem's linear form
				 *  \param residualFormSolutionName a string that identifies the name of the function representing the problem
				 *  solution in the ResidualForm 
				 *  \param jacobianFormSolutionName a string that identifies the name of the function representing the problem
				 *  solution in the JacobianForm. Default value is the empty string, in which case \c residualFormSolutionName_
				 *  will be used for jacobianFormSolutionName_
			 	 *  The stored mesh's and function space's ownership will be shared between the object and the input argument.
			 	 *  The bilinear and linear form will be created too, calling the constructor which takes the function space
			 	 *  as input.
			 	 */
				NonLinearDifferentialProblem (const std::shared_ptr<dolfin::Mesh> mesh, 
			                           	      const std::shared_ptr<dolfin::FunctionSpace> functionSpace,
			                           	      const T_ResidualForm& residualForm,
			                           	      const T_JacobianForm& jacobianForm,
			                           	      const std::string& residualFormSolutionName,
			                           	      const std::string& jacobianFormSolutionName);

				//! Constructor with references [2]
				/*!
			 	 *  \param mesh the problem mesh as a const dolfin::Mesh&
			 	 *  \param functionSpace the problem finite element space as a const dolfin::FunctionSpace&
				 *  \param residualForm a const reference to the problem's bilinear form
				 *  \param jacobianForm a const reference to the problem's linear form
				 *  \param residualFormSolutionName a string that identifies the name of the function representing the problem
				 *  solution in the ResidualForm 
				 *  \param jacobianFormSolutionName a string that identifies the name of the function representing the problem
				 *  solution in the JacobianForm. Default value is the empty string, in which case \c residualFormSolutionName_
				 *  will be used
			 	 *  The stored mesh's and function space's ownership will be unique to the object, since the pointers are 
			 	 *  initialized using the \c new operator and mesh's and functionSpace's copy constructor
			 	 *  The bilinear and linear form will be created too, calling the constructor which takes the function space
			 	 *  as input.
			 	 */
				NonLinearDifferentialProblem (const dolfin::Mesh& mesh, 
			                           	      const dolfin::FunctionSpace& functionSpace,
			                           	      const T_ResidualForm& residualForm,
			                           	      const T_JacobianForm& jacobianForm,
			                           	      const std::string& residualFormSolutionName,
			                           	      const std::string& jacobianFormSolutionName);

				//! Constructor with rvalue references [3]
				/*!
			 	 *  \param mesh the problem mesh as a dolfin::Mesh&&
			 	 *  \param functionSpace the problem finite element space as a dolfin::FunctionSpace&&
				 *  \param residualForm a rvalue reference to the problem's bilinear form
				 *  \param jacobianForm a rvalue reference to the problem's linear form
				 *  \param residualFormSolutionName a string that identifies the name of the function representing the problem
				 *  solution in the ResidualForm 
				 *  \param jacobianFormSolutionName a string that identifies the name of the function representing the problem
				 *  solution in the JacobianForm. Default value is the empty string, in which case \c residualFormSolutionName_
				 *  will be used for jacobianFormSolutionName_
			 	 *  The stored mesh's and function space's ownership will be unique to the object, since the pointers are 
			 	 *  initialized using the \c new operator and mesh's and functionSpace's move constructor
			 	 *  The bilinear and linear form will be created too, calling the constructor which takes the function space
			 	 *  as input.
			 	 */
				NonLinearDifferentialProblem (dolfin::Mesh&& mesh, 
			                           	      dolfin::FunctionSpace&& functionSpace,
			                           	      T_ResidualForm&& residualForm,
			                           	      T_JacobianForm&& jacobianForm,
			                           	      const std::string& residualFormSolutionName,
			                           	      const std::string& jacobianFormSolutionName);

				
				/******************* GETTERS *******************/
				//! Get problem's mesh
				/*! 
			 	 *  \return a const reference to the problem's mesh
			 	 */
				const dolfin::Mesh& mesh () const;

				//! Get problem's finite element space
				/*! 
			 	 *  \return a const reference to the problem's function space
			 	 */
				const dolfin::FunctionSpace& functionSpace () const;

				//! Get const reference to the problem's bilinear form
				/*! 
			 	 *  \return a const reference to the problem's bilinear form
			 	 */
				const T_ResidualForm& residualForm () const;

				//! Get const reference to the problem's linear form
				/*! 
			 	 *  \return a const reference to the problem's linear form
			 	 */
				const T_JacobianForm& jacobianForm () const;

				//! Get const reference to the problem's i-th dirichlet boundary condition
				/*! 
			 	 *  \param i the position in the vector storing all dirichlet BCs of the dirichlet BC object to be 
			 	 *            returned. If called with no argument, \f$i = 0\f$ is assumed. 
			 	 *  \return a const reference to the problem's dirichletBC in position i.
			 	 *          No check is performed on the input value
			 	 */
				const dolfin::DirichletBC& dirichletBC (const std::size_t& i) const;

				//! Get const reference to the problem's dirichlet boundary conditions vector
				/*! 
			 	 *  \return a const reference to the problem's dirichletBC vector
			 	 */
				const std::vector<const dolfin::DirichletBC>& dirichletBC () const;

				//! Get const reference to the solver's parameters
				/*! 
			 	 *  \return a const reference to the solver's parameters
			 	 */
				const dolfin::Parameters& solverParameters () const;

				//! Get const reference to the problem's solution
				/*!
			 	 *  \return a const reference to the problem's solution
			 	 */
				const dolfin::Function& solution () const; 	


				/******************* SETTERS *******************/

				//! Set residual form coefficient [1]
				/*!
			 	 *  Set coefficient with given number. It simply wraps the set_coefficient dolfin function
			 	 *  \param i number of the coefficient to be set
			 	 *  \param coefficient coefficient's value. We use boost::shared_ptr instead of std::shared_ptr for
				 *  compatibility with dolfin1.3.0
			 	 */
				void setResidualFormCoefficient (const std::size_t& i, 
				                                 const boost::shared_ptr<const dolfin::GenericFunction> coefficient);

				//! Set residual form coefficient [2]
				/*!
			 	 *  Set coefficient with given name. It simply wraps the set_coefficient dolfin function
			 	 *  \param name name of the coefficient to be set
			 	 *  \param coefficient coefficient's value. We use boost::shared_ptr instead of std::shared_ptr for
				 *  compatibility with dolfin1.3.0
			 	 */
				void setResidualFormCoefficient (const std::string& name, 
				                                 const boost::shared_ptr<const dolfin::GenericFunction> coefficient);

				//! Set residual form coefficients
				/*!
			 	 *  Set coefficients in given map. It simply wraps the set_coefficient dolfin function
			 	 *  \param coefficients map of the coefficients to be set. 
				 *  We use boost::shared_ptr instead of std::shared_ptr for compatibility with dolfin1.3.0
			 	 */
				void setResidualFormCoefficients 
					(const std::map<std::string, boost::shared_ptr<const dolfin::GenericFunction>>& coefficients);

				//! Set jacobian form coefficient [1]
				/*!
			 	 *  Set coefficient with given number. It simply wraps the set_coefficient dolfin function
			 	 *  \param i number of the coefficient to be set
			 	 *  \param coefficient coefficient's value. We use boost::shared_ptr instead of std::shared_ptr 
				 *  for compatibility with dolfin1.3.0
			 	 */
				void setJacobianFormCoefficient (const std::size_t& i, 
				                                 const boost::shared_ptr<const dolfin::GenericFunction> coefficient);

				//! Set jacobian form coefficient [2]
				/*!
			 	 *  Set coefficient with given name. It simply wraps the set_coefficient dolfin function
			 	 *  \param name name of the coefficient to be set
			 	 *  \param coefficient coefficient's value.
				 *  We use boost::shared_ptr instead of std::shared_ptr for compatibility with dolfin1.3.0
			 	 */
				void setJacobianFormCoefficient (const std::string& name, 
				                                 const boost::shared_ptr<const dolfin::GenericFunction> coefficient);

				//! Set jacobian form coefficients
				/*!
			 	 *  Set coefficients in given map. It simply wraps the set_coefficient dolfin function
			 	 *  \param coefficients map of the coefficients to be set.
				 *  We use boost::shared_ptr instead of std::shared_ptr for compatibility with dolfin1.3.0
			 	 */
				void setJacobianFormCoefficients 
					(const std::map<std::string, boost::shared_ptr<const dolfin::GenericFunction>>& coefficients);
				
				//! Set integration subdomains for the residual form
				/*! Input arguments are:
				 *  \param meshFunction the mesh function used to set the integration subdomains. Template-ized over
				 *  the type of stored objects
				 *  \param subdomainType the type of the subdomains, chosen among those provided by the enumeration
				 *  class \c control_problem::SubdomainType
				 */
				template <class T>
					void setResidualFormIntegrationSubdomains (const dolfin::MeshFunction<T>& meshFunction, 
					                                           const control_problem::SubdomainType& subdomainType);

				//! Set integration subdomains for the jacobian form
				/*! Input arguments are:
				 *  \param meshFunction the mesh function used to set the integration subdomains. Template-ized over
				 *  the type of stored objects
				 *  \param subdomainType the type of the subdomains, chosen among those provided by the enumeration
				 *  class \c SubdomainType
				 */
				template <class T>
					void setJacobianFormIntegrationSubdomains (const dolfin::MeshFunction<T>& meshFunction, 
					                                         const control_problem::SubdomainType& subdomainType);
				
				//! Add Dirichlet boundary condition to the problem [1]
				/*!
			 	 *  \param dirichletCondition a const reference to the dirichlet boundary condition to be added to the problem
			 	 */
				void addDirichletBC (const dolfin::DirichletBC& dirichletCondition);

				//! Add Dirichlet boundary condition to the problem [2]
				/*!
			 	 *  \param dirichletCondition a rvalue reference to the dirichlet boundary condition to be added to the problem
			 	 */
				void addDirichletBC (dolfin::DirichletBC&& dirichletCondition);

				//! Remove Dirichlet boundary condition with given position
				/*!
			 	 *  \param i the position in the vector of the boundary condition to be removed.
			 	 *            If i is greater than the size of the vector, nothing is removed.
			 	 */
				void removeDirichletBC (const std::vector<dolfin::DirichletBC>::iterator& i);

				//! Set solver parameters
				/*!
				 *  \param parameters object of type dolfin::Parameters that contains all the parameters that
				 *  should be passed to the solver
				 */
				void setSolverParameters (const dolfin::Parameters& parameters);

				//! Set initial guess for Newton-Raphson iterative method [1]
				/*!
				 *  \param initialGuess a reference to a dolfin::Function object
				 */
				void setInitialGuess (const dolfin::Function& initialGuess);
				
				//! Set initial guess for Newton-Raphson iterative method [2]
				/*!
				 *  \param initialGuess a reference to a dolfin::Expression object
				 */
				void setInitialGuess (const dolfin::Expression& initialGuess);
				
				
				/******************* METHODS *******************/

				//! Solve problem
				/*!
			 	 * This method solves the problem defined. It uses the private members' value to set the problem and then
			 	 * stores the solution in the private member \c solution_
			 	 */
				virtual void solve ();

				//! Solve problem specifying flag
				/*!
			 	 *  \param solverParameters object of type dolfin::parameters that contain the parameters to
				 *  be used for the non linear solver
			 	 */
				void solve (const dolfin::Parameters& solverParameters);


				/******************* DESTRUCTOR *******************/
				
				//! Destructor
				/*! Default destructor, since members of the class are trivially 
			 	 * destructible.
			 	 * It is declared virtual so that derived classes' constructor
			 	 * can be called on derived classes.
			 	 * The "default-ness" is set in implementation outside of the class for compatibility with
			 	 * gcc-4.6, which does not allow virtual members to be defaulted in class
			 	 */
				virtual ~NonLinearDifferentialProblem ();

				// ---------------------------------------------------------------------------------------------//

			protected:
				//! The problem mesh
				/*! 
			 	 *  Stored as a shared_ptr because it may be common to more than 
			 	 *  one problem
			 	 */
				std::shared_ptr<dolfin::Mesh> mesh_;

				//! The problem finite element space
				/*! 
			 	 *  Stored as a shared_ptr because it may be common to more than 
			 	 *  one problem
			 	 */
				std::shared_ptr<dolfin::FunctionSpace> functionSpace_;

				//! The bilinear form
				T_ResidualForm residualForm_;

				//! The linear form
				T_JacobianForm jacobianForm_;
				
				//! The name of the function representing the solution in ResidualForm
				std::string residualFormSolutionName_;

				//! The name of the function representing the solution in JacobianForm
				std::string jacobianFormSolutionName_;

				//! The Dirichlet's boundary conditions vector
				std::vector<dolfin::DirichletBC> dirichletBCs_;

				//! the solver parameters to be used
				dolfin::Parameters solverParameters_;

				// ---------------------------------------------------------------------------------------------//

			private:
		};



	// ==============================================================================================//
	// ==================================== IMPLEMENTATION ==========================================//
	// ==============================================================================================//


	/******************* CONSTRUCTORS *******************/

	template <class T_ResidualForm, class T_JacobianForm>
		NonLinearDifferentialProblem<T_ResidualForm, T_JacobianForm>::
		NonLinearDifferentialProblem (const std::shared_ptr<dolfin::Mesh> mesh, 
		                              const std::shared_ptr<dolfin::FunctionSpace> functionSpace,
			                          const std::string& residualFormSolutionName,
			                          const std::string& jacobianFormSolutionName = "") : 
			AbstractDifferentialProblem (*functionSpace),
			mesh_ (mesh),
			functionSpace_ (functionSpace),
			residualForm_ (*functionSpace),
			jacobianForm_ (*functionSpace, *functionSpace),
			residualFormSolutionName_ (residualFormSolutionName),
			jacobianFormSolutionName_ (jacobianFormSolutionName.empty () ? residualFormSolutionName : jacobianFormSolutionName),
			dirichletBCs_ (),
			solverParameters_ (dolfin::NonlinearVariationalSolver::default_parameters ())
		{ 
			// set solution in residual and jacobian form
			boost::shared_ptr <dolfin::Function> tmpSolution (new dolfin::Function (solution_));
			residualForm_.set_coefficient (residualFormSolutionName_, tmpSolution);
			jacobianForm_.set_coefficient (jacobianFormSolutionName_, tmpSolution);
		}



	template <class T_ResidualForm, class T_JacobianForm>
		NonLinearDifferentialProblem<T_ResidualForm, T_JacobianForm>::
		NonLinearDifferentialProblem (const dolfin::Mesh& mesh, 
		                              const dolfin::FunctionSpace& functionSpace,
			                          const std::string& residualFormSolutionName,
			                          const std::string& jacobianFormSolutionName = "") : 
			AbstractDifferentialProblem (functionSpace),
			mesh_ (new dolfin::Mesh (mesh)),
			functionSpace_ (new dolfin::FunctionSpace (functionSpace)),
			residualForm_ (functionSpace),
			jacobianForm_ (functionSpace, functionSpace),
			residualFormSolutionName_ (residualFormSolutionName),
			jacobianFormSolutionName_ (jacobianFormSolutionName.empty () ? residualFormSolutionName : jacobianFormSolutionName),
			dirichletBCs_ (),
			solverParameters_ (dolfin::NonlinearVariationalSolver::default_parameters ())
		{ 
			// set solution in residual and jacobian form
			boost::shared_ptr <dolfin::Function> tmpSolution (new dolfin::Function (solution_));
			residualForm_.set_coefficient (residualFormSolutionName_, tmpSolution);
			jacobianForm_.set_coefficient (jacobianFormSolutionName_, tmpSolution);
		}
			

	
	template <class T_ResidualForm, class T_JacobianForm>
		NonLinearDifferentialProblem<T_ResidualForm, T_JacobianForm>::
		NonLinearDifferentialProblem (dolfin::Mesh&& mesh, 
		                              dolfin::FunctionSpace&& functionSpace,
			                          const std::string& residualFormSolutionName,
			                          const std::string& jacobianFormSolutionName = "") : 
			AbstractDifferentialProblem (functionSpace),
			mesh_ (std::make_shared<dolfin::Mesh> (mesh)),
			functionSpace_ (std::make_shared<dolfin::FunctionSpace> (functionSpace)),
			residualForm_ (functionSpace),
			jacobianForm_ (functionSpace, functionSpace),
			residualFormSolutionName_ (residualFormSolutionName),
			jacobianFormSolutionName_ (jacobianFormSolutionName.empty () ? residualFormSolutionName : jacobianFormSolutionName),
			dirichletBCs_ (),
			solverParameters_ (dolfin::NonlinearVariationalSolver::default_parameters ())
		{ 
			// set solution in residual and jacobian form
			boost::shared_ptr <dolfin::Function> tmpSolution (new dolfin::Function (solution_));
			residualForm_.set_coefficient (residualFormSolutionName_, tmpSolution);
			jacobianForm_.set_coefficient (jacobianFormSolutionName_, tmpSolution);
		}
			


	template <class T_ResidualForm, class T_JacobianForm>
		NonLinearDifferentialProblem<T_ResidualForm, T_JacobianForm>::
		NonLinearDifferentialProblem (const std::shared_ptr<dolfin::Mesh> mesh, 
		                              const std::shared_ptr<dolfin::FunctionSpace> functionSpace,
		                              const T_ResidualForm& residualForm,
		                              const T_JacobianForm& jacobianForm,
			                          const std::string& residualFormSolutionName,
			                          const std::string& jacobianFormSolutionName = "") : 
			AbstractDifferentialProblem (*functionSpace),
			mesh_ (mesh),
			functionSpace_ (functionSpace),
			residualForm_ (residualForm),
			jacobianForm_ (jacobianForm),
			residualFormSolutionName_ (residualFormSolutionName),
			jacobianFormSolutionName_ (jacobianFormSolutionName.empty () ? residualFormSolutionName : jacobianFormSolutionName),
			dirichletBCs_ (),
			solverParameters_ (dolfin::NonlinearVariationalSolver::default_parameters ())
		{  }

	

	template <class T_ResidualForm, class T_JacobianForm>
		NonLinearDifferentialProblem<T_ResidualForm, T_JacobianForm>::
		NonLinearDifferentialProblem (const dolfin::Mesh& mesh, 
		                              const dolfin::FunctionSpace& functionSpace,
		                              const T_ResidualForm& residualForm,
		                              const T_JacobianForm& jacobianForm,
			                          const std::string& residualFormSolutionName,
			                          const std::string& jacobianFormSolutionName = "") : 
			AbstractDifferentialProblem (functionSpace),
			mesh_ (new dolfin::Mesh (mesh)),
			functionSpace_ (new dolfin::FunctionSpace (functionSpace)),
			residualForm_ (residualForm),
			jacobianForm_ (jacobianForm),
			residualFormSolutionName_ (residualFormSolutionName),
			jacobianFormSolutionName_ (jacobianFormSolutionName.empty () ? residualFormSolutionName : jacobianFormSolutionName),
			dirichletBCs_ (),
			solverParameters_ (dolfin::NonlinearVariationalSolver::default_parameters ())
		{  }
			


	template <class T_ResidualForm, class T_JacobianForm>
		NonLinearDifferentialProblem<T_ResidualForm, T_JacobianForm>::
		NonLinearDifferentialProblem (dolfin::Mesh&& mesh, 
		                              dolfin::FunctionSpace&& functionSpace,
		                              T_ResidualForm&& residualForm,
		                              T_JacobianForm&& jacobianForm,
			                          const std::string& residualFormSolutionName,
			                          const std::string& jacobianFormSolutionName = "") : 
			AbstractDifferentialProblem (functionSpace),
			mesh_ (std::make_shared<dolfin::Mesh> (mesh)),
			functionSpace_ (std::make_shared<dolfin::FunctionSpace> (functionSpace)),
			residualForm_ (residualForm),
			jacobianForm_ (jacobianForm),
			residualFormSolutionName_ (residualFormSolutionName),
			jacobianFormSolutionName_ (jacobianFormSolutionName.empty () ? residualFormSolutionName : jacobianFormSolutionName),
			dirichletBCs_ (),
			solverParameters_ (dolfin::NonlinearVariationalSolver::default_parameters ())
		{  }
	



	/******************* DESTRUCTOR *******************/

	// this is done for compatibility with gcc-4.6, which doesn't allow virtual members to be defualted in class body
	template <class T_ResidualForm, class T_JacobianForm>
		NonLinearDifferentialProblem<T_ResidualForm, T_JacobianForm>::
		~NonLinearDifferentialProblem () = default;
	
	


	/******************* GETTERS *******************/

	template <class T_ResidualForm, class T_JacobianForm>
		const dolfin::Mesh& NonLinearDifferentialProblem<T_ResidualForm, T_JacobianForm>::
		mesh () const
		{
			return *mesh_;		
		}



	template <class T_ResidualForm, class T_JacobianForm>
		const dolfin::FunctionSpace& NonLinearDifferentialProblem<T_ResidualForm, T_JacobianForm>::
		functionSpace () const
		{
			return *functionSpace_;
		}



	template <class T_ResidualForm, class T_JacobianForm>
		const T_ResidualForm& NonLinearDifferentialProblem<T_ResidualForm, T_JacobianForm>::
		residualForm () const
		{
			return residualForm_;
		}



	template <class T_ResidualForm, class T_JacobianForm>
		const T_JacobianForm& NonLinearDifferentialProblem<T_ResidualForm, T_JacobianForm>::
		jacobianForm () const
		{
			return jacobianForm_;
		}



	template <class T_ResidualForm, class T_JacobianForm>
		const dolfin::DirichletBC& NonLinearDifferentialProblem<T_ResidualForm, T_JacobianForm>::
		dirichletBC (const std::size_t& i) const
		{
			return dirichletBCs_[i];
		}



	template <class T_ResidualForm, class T_JacobianForm>
		const std::vector<const dolfin::DirichletBC>& NonLinearDifferentialProblem<T_ResidualForm, T_JacobianForm>::
		dirichletBC () const
		{
			return dirichletBCs_;
		}



	template <class T_ResidualForm, class T_JacobianForm>
		const dolfin::Parameters& NonLinearDifferentialProblem<T_ResidualForm, T_JacobianForm>::
		solverParameters () const
		{
			return solverParameters_;
		}


	template <class T_ResidualForm, class T_JacobianForm>
		const dolfin::Function& NonLinearDifferentialProblem<T_ResidualForm, T_JacobianForm>::
		solution () const
		{
			return solution_;
		}



	/******************* SETTERS *******************/

	template <class T_ResidualForm, class T_JacobianForm>
		void NonLinearDifferentialProblem<T_ResidualForm, T_JacobianForm>::
		setResidualFormCoefficient (const std::size_t& i, const boost::shared_ptr<const dolfin::GenericFunction> coefficient)
		{
			residualForm_.set_coefficient (i, coefficient);
		}



	template <class T_ResidualForm, class T_JacobianForm>
		void NonLinearDifferentialProblem<T_ResidualForm, T_JacobianForm>::
		setResidualFormCoefficient (const std::string& name, const boost::shared_ptr<const dolfin::GenericFunction> coefficient)
		{
			residualForm_.set_coefficient (name, coefficient);
		}



	template <class T_ResidualForm, class T_JacobianForm>
		void NonLinearDifferentialProblem<T_ResidualForm, T_JacobianForm>::
		setResidualFormCoefficients (const std::map<std::string, boost::shared_ptr<const dolfin::GenericFunction>>& coefficients)
		{
			residualForm_.set_coefficients (coefficients);
		}



	template <class T_ResidualForm, class T_JacobianForm>
		void NonLinearDifferentialProblem<T_ResidualForm, T_JacobianForm>::
		setJacobianFormCoefficient (const std::size_t& i, const boost::shared_ptr<const dolfin::GenericFunction> coefficient)
		{
			jacobianForm_.set_coefficient (i, coefficient);
		}



	template <class T_ResidualForm, class T_JacobianForm>
		void NonLinearDifferentialProblem<T_ResidualForm, T_JacobianForm>::
		setJacobianFormCoefficient (const std::string& name, const boost::shared_ptr<const dolfin::GenericFunction> coefficient)
		{
			jacobianForm_.set_coefficient (name, coefficient);
		}



	template <class T_ResidualForm, class T_JacobianForm>
		void NonLinearDifferentialProblem<T_ResidualForm, T_JacobianForm>::
		setJacobianFormCoefficients (const std::map<std::string, boost::shared_ptr<const dolfin::GenericFunction>>& coefficients)
		{
			jacobianForm_.set_coefficients (coefficients);
		}



	template <class T_ResidualForm, class T_JacobianForm>
		template <class T>
			void NonLinearDifferentialProblem<T_ResidualForm, T_JacobianForm>::
			setResidualFormIntegrationSubdomains (const dolfin::MeshFunction<T>& meshFunction, 
			                                      const control_problem::SubdomainType& subdomainType)
			{
				if (subdomainType == control_problem::SubdomainType::INTERNAL_CELLS)
				{
					residualForm_.dx = meshFunction;
				}
				else if (subdomainType == control_problem::SubdomainType::INTERNAL_FACETS)
				{
					residualForm_.dS = meshFunction;
				}
				else if (subdomainType == control_problem::SubdomainType::BOUNDARY_FACETS)
				{
					residualForm_.ds = meshFunction;
				}
				else
				{
					dolfin::warning ("unknown subdomain type requested while trying to apply mesh function to residual form");
				}
			}



	template <class T_ResidualForm, class T_JacobianForm>
		template <class T>
			void NonLinearDifferentialProblem<T_ResidualForm, T_JacobianForm>::
			setJacobianFormIntegrationSubdomains (const dolfin::MeshFunction<T>& meshFunction, 
		                                    	const control_problem::SubdomainType& subdomainType)
			{
				if (subdomainType == control_problem::SubdomainType::INTERNAL_CELLS)
					{
						jacobianForm_.dx = meshFunction;
					}
					else if (subdomainType == control_problem::SubdomainType::INTERNAL_FACETS)
					{
						jacobianForm_.dS = meshFunction;
					}
					else if (subdomainType == control_problem::SubdomainType::BOUNDARY_FACETS)
					{
						jacobianForm_.ds = meshFunction;
					}
					else
					{
						dolfin::warning ("unknown subdomain type requested while trying to apply mesh function to jacobian form");
					}
				}



	template <class T_ResidualForm, class T_JacobianForm>
		void NonLinearDifferentialProblem<T_ResidualForm, T_JacobianForm>::
		addDirichletBC (const dolfin::DirichletBC& dirichletCondition)
		{
			dirichletBCs_.emplace_back (dirichletCondition);
		}



	template <class T_ResidualForm, class T_JacobianForm>
		void NonLinearDifferentialProblem<T_ResidualForm, T_JacobianForm>::
		addDirichletBC (dolfin::DirichletBC&& dirichletCondition)
		{
			dirichletBCs_.emplace_back (dirichletCondition);
		}



	template <class T_ResidualForm, class T_JacobianForm>
		void NonLinearDifferentialProblem<T_ResidualForm, T_JacobianForm>::
		removeDirichletBC (const std::vector<dolfin::DirichletBC>::iterator& i)
		{
			dirichletBCs_.erase (i);
		}



	template <class T_ResidualForm, class T_JacobianForm>
		void NonLinearDifferentialProblem<T_ResidualForm, T_JacobianForm>::
		setSolverParameters (const dolfin::Parameters& parameters)
		{
			solverParameters_.update (parameters);
		}



	template <class T_ResidualForm, class T_JacobianForm>
		void NonLinearDifferentialProblem<T_ResidualForm, T_JacobianForm>::
		setInitialGuess (const dolfin::Function& initialGuess)
		{
			solution_ = initialGuess;
		}



	template <class T_ResidualForm, class T_JacobianForm>
		void NonLinearDifferentialProblem<T_ResidualForm, T_JacobianForm>::
		setInitialGuess (const dolfin::Expression& initialGuess)
		{
			solution_ = initialGuess;
		}



	template <class T_ResidualForm, class T_JacobianForm>
		void NonLinearDifferentialProblem<T_ResidualForm, T_JacobianForm>::
		solve () 
		{
			// create vector of pointers to DirichletBC. This is needed to call the function dolfin::solve
			std::vector<const dolfin::DirichletBC*> tmpDirichletBCs (dirichletBCs_.size (), nullptr);
			for (std::size_t i = 0; i < dirichletBCs_.size (); ++i)
			{
				tmpDirichletBCs[i] = &dirichletBCs_[i];
			}
			// note that when tmpDirichletBCs gets destroyd (upon exit of the function) the vector distructor will
			// call DirichletBC* destructor, which means that the pointers will be destroyed but the object they point
			// to will not
			
			// now solve non linear problem and store solution in solution_
			// solverParameters_ might be empty. It can be filled using the appropriate class functions
			// we need to check whether tmpDirichletBCs is empty, and call the appropriate function
			if (tmpDirichletBCs.size () != 0)
			{
				dolfin::solve (residualForm_ == 0, solution_, tmpDirichletBCs, jacobianForm_, solverParameters_);
			}
			else
			{
				dolfin::solve (residualForm_ == 0, solution_, jacobianForm_, solverParameters_);
			}
		}



	template <class T_ResidualForm, class T_JacobianForm>
		void NonLinearDifferentialProblem<T_ResidualForm, T_JacobianForm>::
		solve (const dolfin::Parameters& solverParameters)
		{
			// create vector of pointers to DirichletBC. This is needed to call the function dolfin::solve
			std::vector<const dolfin::DirichletBC*> tmpDirichletBCs (dirichletBCs_.size (), nullptr);
			for (std::size_t i = 0; i < dirichletBCs_.size (); ++i)
			{
				tmpDirichletBCs[i] = &dirichletBCs_[i];
			}
			// note that when tmpDirichletBCs gets destroyd (upon exit of the function) the vector distructor will
			// call DirichletBC* destructor, which means that the pointers will be destroyed but the object they point
			// to will not
			
			// now solve non linear problem and store solution in solution_
			// solverParameters_ might be empty. It can be filled using the appropriate class functions
			// we need to check whether tmpDirichletBCs is empty, and call the appropriate function
			if (tmpDirichletBCs.size () != 0)
			{
				dolfin::solve (residualForm_ == 0, solution_, tmpDirichletBCs, jacobianForm_, solverParameters);
			}
			else
			{
				dolfin::solve (residualForm_ == 0, solution_, jacobianForm_, solverParameters);
			}
		}
}
#endif

