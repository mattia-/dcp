#ifndef SRC_DIFFERENTIALPROBLEM_LINEARDIFFERENTIALPROBLEM_HPP_INCLUDE_GUARD
#define SRC_DIFFERENTIALPROBLEM_LINEARDIFFERENTIALPROBLEM_HPP_INCLUDE_GUARD

#include <dolfin.h>
#include <vector>
#include <map>
#include <string>
#include <memory>
#include <DifferentialProblem/AbstractDifferentialProblem.hpp>
#include <Factory/LinearSolverFactory.hpp>
#include <Utils/SubdomainType.hpp>

namespace control_problem
{
	/*! \class LinearDifferentialProblem LinearDifferentialProblem.hpp
	 *  \brief Class for linear differential problems.
	 *
	 *  This class represents problem of the form
	 *  \f[
	 *      \mbox{Find } u \in V : a \left(u, v\right) = F \left(v\right) \ \forall\,v\,\in\,V
	 *  \f]
	 *  with \f$ a \left(u, v\right) : V \times V \rightarrow \mathds{R}\f$ bilinear form on \f$V\f$
	 *  and \f$ L \left(v\right) : V \rightarrow \mathds{R} \f$ linear form on the same space.
	 *  
	 *  It inherits publicly from \c AbstractDifferentialProblem
	 *  and it extends its functionalities to a concrete differential
	 *  problem.
	 *  Template arguments are:
	 *  \arg T_BilinearForm the bilinear form type
	 *  \arg T_LinearForm the linear form type
	 *  \arg T_LinearSolverFactory the type of the factory that creates the linear solver. By default, it is set
	 *  to control_problem::LinearSolverFactory
	 *  This is needed because the dolfin function "coefficient_number" was neither declared
	 *  virtual nor implemented by FEniCS's developers in the dolfin class "Form", from which
	 *  we first thought to derive for the protected members of the LinearDifferentialProblem
	 *  class
	 */

	template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory = control_problem::LinearSolverFactory>
		class LinearDifferentialProblem : public AbstractDifferentialProblem
		{
			// ---------------------------------------------------------------------------------------------//	

			public:
				/******************* CONSTRUCTORS *******************/
				//! Default constructor is deleted. The class is not default constructable.
				LinearDifferentialProblem () = delete;

				//!  Constructor with shared pointers [1]
				/*!
			 	 *  \param mesh the problem mesh as a const std::shared_ptr to dolfin::Mesh
			 	 *  \param functionSpace the problem finite element space as a const std::shared_ptr to dolfin::FunctionSpace
				 *  \param solverType the type of the solver. Default: lu_solver
				 *  \param solverMethod the method of the solver. Possible values depend on the solver type (see dolfin
				 *  documentation, or use method \c list_<solverType>_methods). Default value: default
				 *  \param solverPreconditioner the preconditioner to be used. It is not used for lu_solvers. Default
				 *  value: default
			 	 *  The stored mesh's and function space's ownership will be shared between the object and the input argument.
			 	 *  The bilinear and linear form will be created too, calling the constructor which takes the function space
			 	 *  as input.
			 	 */
				LinearDifferentialProblem (const std::shared_ptr<dolfin::Mesh> mesh, 
			                           	   const std::shared_ptr<dolfin::FunctionSpace> functionSpace,
			                           	   const std::string& solverType,
			                           	   const std::string& solverMethod,
			                           	   const std::string& solverPreconditioner);
				

				//! Constructor with references [1]
				/*!
			 	 *  \param mesh the problem mesh as a const dolfin::Mesh&
			 	 *  \param functionSpace the problem finite element space as a const dolfin::FunctionSpace&
				 *  \param solverType the type of the solver. Default: lu_solver
				 *  \param solverMethod the method of the solver. Possible values depend on the solver type (see dolfin
				 *  documentation, or use method \c list_<solverType>_methods). Default value: default
				 *  \param solverPreconditioner the preconditioner to be used. It is not used for lu_solvers. Default
				 *  value: default
			 	 *  The stored mesh's and function space's ownership will be unique to the object, since the pointers are 
			 	 *  initialized using the \c new operator and mesh's and functionSpace's copy constructor
			 	 *  The bilinear and linear form will be created too, calling the constructor which takes the function space
			 	 *  as input.
			 	 */
				LinearDifferentialProblem (const dolfin::Mesh& mesh, 
			                           	   const dolfin::FunctionSpace& functionSpace,
			                           	   const std::string& solverType,
			                           	   const std::string& solverMethod,
			                           	   const std::string& solverPreconditioner);

				//! Constructor with rvalue references [1]
				/*!
			 	 *  \param mesh the problem mesh as a dolfin::Mesh&&
			 	 *  \param functionSpace the problem finite element space as a dolfin::FunctionSpace&&
				 *  \param solverType the type of the solver. Default: lu_solver
				 *  \param solverMethod the method of the solver. Possible values depend on the solver type (see dolfin
				 *  documentation, or use method \c list_<solverType>_methods). Default value: default
				 *  \param solverPreconditioner the preconditioner to be used. It is not used for lu_solvers. Default
				 *  value: default
			 	 *  The stored mesh's and function space's ownership will be unique to the object, since the pointers are 
			 	 *  initialized using the \c new operator and mesh's and functionSpace's move constructor
			 	 *  The bilinear and linear form will be created too, calling the constructor which takes the function space
			 	 *  as input.
			 	 */
				LinearDifferentialProblem (dolfin::Mesh&& mesh, 
			                           	   dolfin::FunctionSpace&& functionSpace,
			                           	   const std::string& solverType,
			                           	   const std::string& solverMethod,
			                           	   const std::string& solverPreconditioner);

				
				//!  Constructor with shared pointers [2]
				/*!
			 	 *  \param mesh the problem mesh as a const std::shared_ptr to dolfin::Mesh
			 	 *  \param functionSpace the problem finite element space as a const std::shared_ptr to dolfin::FunctionSpace
				 *  \param bilinearForm a const reference to the problem's bilinear form
				 *  \param linearForm a const reference to the problem's linear form
				 *  \param solverType the type of the solver. Default: lu_solver
				 *  \param solverMethod the method of the solver. Possible values depend on the solver type (see dolfin
				 *  documentation, or use method \c list_<solverType>_methods. Default value: default
				 *  \param solverPreconditioner the preconditioner to be used. It is not used for lu_solvers. Default
				 *  value: default
			 	 *  The stored mesh's and function space's ownership will be shared between the object and the input argument.
			 	 *  The bilinear and linear form will be created too, calling the constructor which takes the function space
			 	 *  as input.
			 	 */
				LinearDifferentialProblem (const std::shared_ptr<dolfin::Mesh> mesh, 
			                           	   const std::shared_ptr<dolfin::FunctionSpace> functionSpace,
			                           	   const T_BilinearForm& bilinearForm,
			                           	   const T_LinearForm& linearForm,
			                           	   const std::string& solverType,
			                           	   const std::string& solverMethod,
			                           	   const std::string& solverPreconditioner);

				//! Constructor with references [2]
				/*!
			 	 *  \param mesh the problem mesh as a const dolfin::Mesh&
			 	 *  \param functionSpace the problem finite element space as a const dolfin::FunctionSpace&
				 *  \param bilinearForm a const reference to the problem's bilinear form
				 *  \param linearForm a const reference to the problem's linear form
				 *  \param solverType the type of the solver. Default: lu_solver
				 *  \param solverMethod the method of the solver. Possible values depend on the solver type (see dolfin
				 *  documentation, or use method \c list_<solverType>_methods. Default value: default
				 *  \param solverPreconditioner the preconditioner to be used. It is not used for lu_solvers. Default
				 *  value: default
			 	 *  The stored mesh's and function space's ownership will be unique to the object, since the pointers are 
			 	 *  initialized using the \c new operator and mesh's and functionSpace's copy constructor
			 	 *  The bilinear and linear form will be created too, calling the constructor which takes the function space
			 	 *  as input.
			 	 */
				LinearDifferentialProblem (const dolfin::Mesh& mesh, 
			                           	   const dolfin::FunctionSpace& functionSpace,
			                           	   const T_BilinearForm& bilinearForm,
			                           	   const T_LinearForm& linearForm,
			                           	   const std::string& solverType,
			                           	   const std::string& solverMethod,
			                           	   const std::string& solverPreconditioner);

				//! Constructor with rvalue references [3]
				/*!
			 	 *  \param mesh the problem mesh as a dolfin::Mesh&&
			 	 *  \param functionSpace the problem finite element space as a dolfin::FunctionSpace&&
				 *  \param bilinearForm a rvalue reference to the problem's bilinear form
				 *  \param linearForm a rvalue reference to the problem's linear form
				 *  \param solverType the type of the solver. Default: lu_solver
				 *  \param solverMethod the method of the solver. Possible values depend on the solver type (see dolfin
				 *  documentation, or use method \c list_<solverType>_methods. Default value: default
				 *  \param solverPreconditioner the preconditioner to be used. It is not used for lu_solvers. Default
				 *  value: default
			 	 *  The stored mesh's and function space's ownership will be unique to the object, since the pointers are 
			 	 *  initialized using the \c new operator and mesh's and functionSpace's move constructor
			 	 *  The bilinear and linear form will be created too, calling the constructor which takes the function space
			 	 *  as input.
			 	 */
				LinearDifferentialProblem (dolfin::Mesh&& mesh, 
			                           	   dolfin::FunctionSpace&& functionSpace,
			                           	   T_BilinearForm&& bilinearForm,
			                           	   T_LinearForm&& linearForm,
			                           	   const std::string& solverType,
			                           	   const std::string& solverMethod,
			                           	   const std::string& solverPreconditioner);

				
				/******************* GETTERS *******************/
				//! Get const reference to the problem's linear form
				/*! 
			 	 *  \return a const reference to the problem's linear form
			 	 */
				const T_BilinearForm& bilinearForm () const;

				//! Get const reference to the problem's linear form
				/*! 
			 	 *  \return a const reference to the problem's linear form
			 	 */
				const T_LinearForm& linearForm () const;

				//! Get const reference to the problem's solver 
				/*! 
			 	 *  \return a const reference to the problem's solver.
			 	 */
				const dolfin::GenericLinearSolver& solver () const;

				//! Get const reference to the problem's linear operator
				/*!
			 	 *  \return a const reference to the problem's linear operator, which
			 	 *  is a dolfin::Matrix
			 	 */
				const dolfin::Matrix& linearOperator () const;

				//! Get const reference to the problem's right hand side
				/*!
			 	 *  \return a const reference to the problem's right hand side, which
			 	 *  is a dolfin::Vector
			 	 */
				const dolfin::Vector& rhs () const;


				/******************* SETTERS *******************/
				
				//! Set coefficient [1]. Override of virtual function in \c AbstractDifferentialProblem.
				/*!
				 *  Possible values for \c coefficientType are:
				 *  \li bilinear_form to set the coefficient in the bilinear form
				 *  \li linear_form to set the coefficient in the linear form
				 *  
				 *  See \c AbstractDifferentialProblem documentation for more details on the function
				 */
				virtual void setCoefficient (const std::string& coefficientType, 
				                             const boost::shared_ptr<const dolfin::GenericFunction> coefficientValue,
				                             const std::string& coefficientName);

				//! Set coefficient [2]. Override of virtual function in \c AbstractDifferentialProblem.
				/*!
				 *  Possible values for \c coefficientType are:
				 *  \li bilinear_form to set the coefficient in the bilinear form
				 *  \li linear_form to set the coefficient in the linear form
				 *  
				 *  See \c AbstractDifferentialProblem documentation for more details on the function
				 */
				virtual void setCoefficient (const std::string& coefficientType,
				                             const boost::shared_ptr<const dolfin::GenericFunction> coefficientValue,
				                             const std::size_t& coefficientNumber);

				//! Set integration subdomains for the forms. Override of virtual function in \c AbstractDifferentialProblem
				/*! 
				 *  Possible values for \c coefficientType are:
				 *  \li bilinear_form to set the integration subdomain in the bilinear form
				 *  \li linear_form to set the integration subdomain in the linear form
				 *  
				 *  See \c AbstractDifferentialProblem documentation for more details on the function
				 */
				virtual	void setIntegrationSubdomains (const std::string& coefficientType,
													   const dolfin::MeshFunction<std::size_t>& meshFunction,
													   const control_problem::SubdomainType& subdomainType);

				
				/******************* METHODS *******************/

				//! Solve problem
				/*!
			 	 *  This method solves the problem defined. It uses the private members' value to set the problem and then
			 	 *  stores the solution in the private member \c solution_. Note that it checks the protected member
				 *  \c parameters_ to decided whether problem's matrix and vector should be reassembled and if
				 *  the values of the parameters "desired_solver_type", "desired_solver_method" and 
				 *  "desired_solver_preconditioner" match the values of "current_solver_type", "current_solver_method"
				 *  and "current_solver_preconditioner". If they differ, it calls \c createSolver()
			 	 */
				virtual void solve ();

				//! Solve problem specifying flag
				/*!
			 	 *  \param mustReassemble true if the system operators (matrix and right hand side vector)
			 	 *         should be reassembled. It is false by default.
			 	 */
				void solve (const bool& mustReassemble);


				/******************* DESTRUCTOR *******************/
				
				//! Destructor
				/*! Default destructor, since members of the class are trivially 
			 	 * destructible.
			 	 * It is declared virtual so that derived classes' constructor
			 	 * can be called on derived classes.
			 	 * The "default-ness" is set in implementation outside of the class for compatibility with
			 	 * gcc-4.6, which does not allow virtual members to be defaulted in class
			 	 */
				virtual ~LinearDifferentialProblem ();

				// ---------------------------------------------------------------------------------------------//

			protected:
				//! Creates a linear solver object of type passed as input. 
				/*!
				 *  The solver will be created using the parameters "desired_solver_type", "desired_solver_method"
				 *  and "desired_solver_preconditioner" set in the protected member \c parameters_. It will
				 *  also set the parameters "current_solver_type", "current_solver_method" and 
				 *  "current_solver_preconditioner" in the same set of parameters substituting the current values with 
				 *  the values being used to create the solver
				 */
				std::unique_ptr<dolfin::GenericLinearSolver> createSolver ();
				
				//! The bilinear form
				T_BilinearForm bilinearForm_;

				//! The linear form
				T_LinearForm linearForm_;

				//! The solver
				/*! 
			 	 *  We use a pointer so that polymorphism can be applied.
			 	 */
				std::unique_ptr<dolfin::GenericLinearSolver> solver_;
				
				//! Matrix to hold the problem's discrete operator
				/*!
			 	 *  We use a boost::shared_ptr for compatibility with FEniCS, version 1.3.0.
			 	 *  Note that there is no write access to this variable from outside the class,
			 	 *  so it is guaranteed to be a unique shared_ptr
			 	 */
				boost::shared_ptr<dolfin::Matrix> problemMatrix_;

				//! Vector to store the right hand side of the discrete problem
				dolfin::Vector rhsVector_;

				// ---------------------------------------------------------------------------------------------//

			private:
		};



	// ============================================================================================== //
	// ==================================== IMPLEMENTATION ========================================== //
	// ============================================================================================== //


	/******************* CONSTRUCTORS *******************/

	template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
		LinearDifferentialProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
		LinearDifferentialProblem (const std::shared_ptr<dolfin::Mesh> mesh, 
		                           const std::shared_ptr<dolfin::FunctionSpace> functionSpace,
			                       const std::string& solverType = "lu_solver",
			                       const std::string& solverMethod = "default",
			                       const std::string& solverPreconditioner = "default") :
			AbstractDifferentialProblem (mesh, functionSpace),
			bilinearForm_ (*functionSpace, *functionSpace),
			linearForm_ (*functionSpace),
			solver_ (nullptr),
			problemMatrix_ (new dolfin::Matrix),
			rhsVector_ ()
		{ 
			dolfin::log (dolfin::DBG, "Building LinearDifferentialProblem...");
			
			dolfin::log (dolfin::DBG, "Setting up parameters...");
			parameters_.add ("current_solver_type", solverType);
			parameters_.add ("current_solver_method", solverMethod);
			parameters_.add ("current_solver_preconditioner", solverPreconditioner);
			parameters_.add ("desired_solver_type", solverType);
			parameters_.add ("desired_solver_method", solverMethod);
			parameters_.add ("desired_solver_preconditioner", solverPreconditioner);
			
			dolfin::log (dolfin::DBG, "Creating solver...");
			solver_ = createSolver ();
			parameters_.add (solver_ -> parameters);
			parameters_.add ("system_is_assembled", false);
			parameters_.add ("force_reassemble_system", false);	
			
			dolfin::log (dolfin::DBG, "LinearDifferentialProblem created!");
		}



	template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
		LinearDifferentialProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
		LinearDifferentialProblem (const dolfin::Mesh& mesh, 
		                           const dolfin::FunctionSpace& functionSpace,
			                       const std::string& solverType = "lu_solver",
			                       const std::string& solverMethod = "default",
			                       const std::string& solverPreconditioner = "default") :
			AbstractDifferentialProblem (mesh, functionSpace),
			bilinearForm_ (functionSpace, functionSpace),
			linearForm_ (functionSpace),
			solver_ (nullptr),
			problemMatrix_ (new dolfin::Matrix),
			rhsVector_ ()
		{ 
			dolfin::log (dolfin::DBG, "Building LinearDifferentialProblem...");
			
			dolfin::log (dolfin::DBG, "Setting up parameters...");
			parameters_.add ("current_solver_type", solverType);
			parameters_.add ("current_solver_method", solverMethod);
			parameters_.add ("current_solver_preconditioner", solverPreconditioner);
			parameters_.add ("desired_solver_type", solverType);
			parameters_.add ("desired_solver_method", solverMethod);
			parameters_.add ("desired_solver_preconditioner", solverPreconditioner);
			
			dolfin::log (dolfin::DBG, "Creating solver...");
			solver_ = createSolver ();
			parameters_.add (solver_ -> parameters);
			parameters_.add ("system_is_assembled", false);
			parameters_.add ("force_reassemble_system", false);	
			
			dolfin::log (dolfin::DBG, "LinearDifferentialProblem created!");
		}



	template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
		LinearDifferentialProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
		LinearDifferentialProblem (dolfin::Mesh&& mesh, 
		                           dolfin::FunctionSpace&& functionSpace,
			                       const std::string& solverType = "lu_solver",
			                       const std::string& solverMethod = "default",
			                       const std::string& solverPreconditioner = "default") :
			AbstractDifferentialProblem (mesh, functionSpace),
			bilinearForm_ (functionSpace, functionSpace),
			linearForm_ (functionSpace),
			solver_ (nullptr),
			problemMatrix_ (new dolfin::Matrix),
			rhsVector_ ()
		{ 
			dolfin::log (dolfin::DBG, "Building LinearDifferentialProblem...");
			
			dolfin::log (dolfin::DBG, "Setting up parameters...");
			parameters_.add ("current_solver_type", solverType);
			parameters_.add ("current_solver_method", solverMethod);
			parameters_.add ("current_solver_preconditioner", solverPreconditioner);
			parameters_.add ("desired_solver_type", solverType);
			parameters_.add ("desired_solver_method", solverMethod);
			parameters_.add ("desired_solver_preconditioner", solverPreconditioner);
		
			dolfin::log (dolfin::DBG, "Creating solver...");
			solver_ = createSolver ();
			parameters_.add (solver_ -> parameters);
			parameters_.add ("system_is_assembled", false);
			parameters_.add ("force_reassemble_system", false);
			
			dolfin::log (dolfin::DBG, "LinearDifferentialProblem created!");
		}



	template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
		LinearDifferentialProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
		LinearDifferentialProblem (const std::shared_ptr<dolfin::Mesh> mesh, 
		                           const std::shared_ptr<dolfin::FunctionSpace> functionSpace,
		                           const T_BilinearForm& bilinearForm,
		                           const T_LinearForm& linearForm,
		                           const std::string& solverType = "lu_solver",
		                           const std::string& solverMethod = "default",
		                           const std::string& solverPreconditioner = "default") :
			AbstractDifferentialProblem (mesh, functionSpace),
			bilinearForm_ (bilinearForm),
			linearForm_ (linearForm),
			solver_ (nullptr),
			problemMatrix_ (new dolfin::Matrix),
			rhsVector_ ()
		{ 
			dolfin::log (dolfin::DBG, "Building LinearDifferentialProblem...");
			
			dolfin::log (dolfin::DBG, "Setting up parameters...");
			parameters_.add ("current_solver_type", solverType);
			parameters_.add ("current_solver_method", solverMethod);
			parameters_.add ("current_solver_preconditioner", solverPreconditioner);
			parameters_.add ("desired_solver_type", solverType);
			parameters_.add ("desired_solver_method", solverMethod);
			parameters_.add ("desired_solver_preconditioner", solverPreconditioner);
			
			dolfin::log (dolfin::DBG, "Creating solver...");
			solver_ = createSolver ();
			parameters_.add (solver_ -> parameters);
			parameters_.add ("system_is_assembled", false);
			parameters_.add ("force_reassemble_system", false);
			
			dolfin::log (dolfin::DBG, "LinearDifferentialProblem created!");
		}



	template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
		LinearDifferentialProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
		LinearDifferentialProblem (const dolfin::Mesh& mesh, 
		                           const dolfin::FunctionSpace& functionSpace,
		                           const T_BilinearForm& bilinearForm,
		                           const T_LinearForm& linearForm,
			                       const std::string& solverType = "lu_solver",
			                       const std::string& solverMethod = "default",
			                       const std::string& solverPreconditioner = "default") :
			AbstractDifferentialProblem (mesh, functionSpace),
			bilinearForm_ (bilinearForm),
			linearForm_ (linearForm),
			solver_ (nullptr),
			problemMatrix_ (new dolfin::Matrix),
			rhsVector_ ()
		{ 
			dolfin::log (dolfin::DBG, "Building LinearDifferentialProblem...");
			
			dolfin::log (dolfin::DBG, "Setting up parameters...");
			parameters_.add ("current_solver_type", solverType);
			parameters_.add ("current_solver_method", solverMethod);
			parameters_.add ("current_solver_preconditioner", solverPreconditioner);
			parameters_.add ("desired_solver_type", solverType);
			parameters_.add ("desired_solver_method", solverMethod);
			parameters_.add ("desired_solver_preconditioner", solverPreconditioner);
			
			dolfin::log (dolfin::DBG, "Creating solver...");
			solver_ = createSolver ();
			parameters_.add (solver_ -> parameters);
			parameters_.add ("system_is_assembled", false);
			parameters_.add ("force_reassemble_system", false);	
			
			dolfin::log (dolfin::DBG, "LinearDifferentialProblem created!");
		}



	template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
		LinearDifferentialProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
		LinearDifferentialProblem (dolfin::Mesh&& mesh, 
		                           dolfin::FunctionSpace&& functionSpace,
		                           T_BilinearForm&& bilinearForm,
		                           T_LinearForm&& linearForm,
		                           const std::string& solverType = "lu_solver",
		                           const std::string& solverMethod = "default",
		                           const std::string& solverPreconditioner = "default") :
			AbstractDifferentialProblem (mesh, functionSpace),
			bilinearForm_ (std::move (bilinearForm)),
			linearForm_ (std::move (linearForm)),
			solver_ (nullptr),
			problemMatrix_ (new dolfin::Matrix),
			rhsVector_ ()
		{
			dolfin::log (dolfin::DBG, "Building LinearDifferentialProblem...");
			
			dolfin::log (dolfin::DBG, "Setting up parameters...");
			parameters_.add ("current_solver_type", solverType);
			parameters_.add ("current_solver_method", solverMethod);
			parameters_.add ("current_solver_preconditioner", solverPreconditioner);
			parameters_.add ("desired_solver_type", solverType);
			parameters_.add ("desired_solver_method", solverMethod);
			parameters_.add ("desired_solver_preconditioner", solverPreconditioner);
			
			dolfin::log (dolfin::DBG, "Creating solver...");
			solver_ = createSolver ();
			parameters_.add (solver_ -> parameters);
			parameters_.add ("system_is_assembled", false);
			parameters_.add ("force_reassemble_system", false);	
			
			dolfin::log (dolfin::DBG, "LinearDifferentialProblem created!");
		}

	
	/******************* DESTRUCTOR *******************/

	// this is done for compatibility with gcc-4.6, which doesn't allow virtual members to be defualted in class body
	template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
		LinearDifferentialProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
		~LinearDifferentialProblem () = default;
	
	


	/******************* GETTERS *******************/

	template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
		const T_BilinearForm& LinearDifferentialProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
		bilinearForm () const
		{
			return bilinearForm_;
		}



	template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
		const T_LinearForm& LinearDifferentialProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
		linearForm () const
		{
			return linearForm_;
		}



	template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
		const dolfin::GenericLinearSolver& LinearDifferentialProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
		solver () const
		{
			return *solver_;
		}



	template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
		const dolfin::Matrix& LinearDifferentialProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
		linearOperator () const
		{
			return *problemMatrix_;
		}



	template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
		const dolfin::Vector& LinearDifferentialProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
		rhs () const
		{
			return rhsVector_;
		}



	/******************* SETTERS *******************/

	template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
		void LinearDifferentialProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
		setCoefficient (const std::string& coefficientType, 
		                const boost::shared_ptr<const dolfin::GenericFunction> coefficientValue,
		                const std::string& coefficientName)
		{
			if (coefficientType == "bilinear_form")
			{
				dolfin::log (dolfin::DBG, "Setting bilinear form coefficient %s...", coefficientName.c_str ());
				bilinearForm_.set_coefficient (coefficientName, coefficientValue);
				dolfin::log (dolfin::DBG, "done!");
			}
			else if (coefficientType == "linear_form")
			{
				dolfin::log (dolfin::DBG, "Setting linear form coefficient %s...", coefficientName.c_str ());
				linearForm_.set_coefficient (coefficientName, coefficientValue);
				dolfin::log (dolfin::DBG, "done!");
			}
			else
			{
				dolfin::warning ("Cannot set coefficient in linear differential problem. Form type %s unknown",
				                 coefficientType.c_str ());
			}

		}

	

	template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
		void LinearDifferentialProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
		setCoefficient (const std::string& coefficientType,
		                const boost::shared_ptr<const dolfin::GenericFunction> coefficientValue,
		                const std::size_t& coefficientNumber)
		{
			if (coefficientType == "bilinear_form")
			{
				dolfin::log (dolfin::DBG, "Setting bilinear form coefficient number %d...", coefficientNumber);
				bilinearForm_.set_coefficient (coefficientNumber, coefficientValue);
				dolfin::log (dolfin::DBG, "done!");
			}
			else if (coefficientType == "linear_form")
			{
				dolfin::log (dolfin::DBG, "Setting linear form coefficient number %d...", coefficientNumber);
				linearForm_.set_coefficient (coefficientNumber, coefficientValue);
				dolfin::log (dolfin::DBG, "done!");
			}
			else
			{
				dolfin::warning ("Cannot set coefficient in linear differential problem. Form type %s unknown",
				                 coefficientType.c_str ());
			}
		}



	template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
		void LinearDifferentialProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
		setIntegrationSubdomains (const std::string& coefficientType,
		                          const dolfin::MeshFunction<std::size_t>& meshFunction,
		                          const control_problem::SubdomainType& subdomainType)
		{
			if (coefficientType == "bilinear_form")
			{
				if (subdomainType == control_problem::SubdomainType::INTERNAL_CELLS)
				{
					dolfin::log (dolfin::DBG, "Setting bilinear form integration subdomain on INTERNAL_CELLS...");
					bilinearForm_.dx = meshFunction;
					dolfin::log (dolfin::DBG, "done!");
				}
				else if (subdomainType == control_problem::SubdomainType::INTERNAL_FACETS)
				{
					dolfin::log (dolfin::DBG, "Setting bilinear form integration subdomain on INTERNAL_FACETS...");
					bilinearForm_.dS = meshFunction;
					dolfin::log (dolfin::DBG, "done!");
				}
				else if (subdomainType == control_problem::SubdomainType::BOUNDARY_FACETS)
				{
					dolfin::log (dolfin::DBG, "Setting bilinear form integration subdomain on BOUNDARY_FACETS...");
					bilinearForm_.ds = meshFunction;
					dolfin::log (dolfin::DBG, "done!");
				}
				else
				{
					dolfin::warning ("unknown subdomain type requested while trying to apply mesh function to bilinear form"); 
				}
			}
			else if (coefficientType == "linear_form")
			{
				if (subdomainType == control_problem::SubdomainType::INTERNAL_CELLS)
				{
					dolfin::log (dolfin::DBG, "Setting linear form integration subdomain on INTERNAL_CELLS...");
					linearForm_.dx = meshFunction;
					dolfin::log (dolfin::DBG, "done!");
				}
				else if (subdomainType == control_problem::SubdomainType::INTERNAL_FACETS)
				{
					dolfin::log (dolfin::DBG, "Setting linear form integration subdomain on INTERNAL_FACETS...");
					linearForm_.dS = meshFunction;
					dolfin::log (dolfin::DBG, "done!");
				}
				else if (subdomainType == control_problem::SubdomainType::BOUNDARY_FACETS)
				{
					dolfin::log (dolfin::DBG, "Setting linear form integration subdomain on BOUNDARY_FACETS...");
					linearForm_.ds = meshFunction;
					dolfin::log (dolfin::DBG, "done!");
				}
				else
				{
					dolfin::warning ("unknown subdomain type requested while trying to apply mesh function to linear form"); 
				}
			}
			else
			{
				dolfin::warning ("Cannot set integration subdomain in linear differential problem. Form type %s unknown",
				                 coefficientType.c_str ());
			}

		}



	template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
		void LinearDifferentialProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
		solve () 
		{
			// define auxiliary string variables
			std::string desiredSolverType = parameters_ ["desired_solver_type"];
			std::string desiredSolverMethod = parameters_ ["desired_solver_method"];
			std::string desiredSolverPreconditioner = parameters_ ["desired_solver_preconditioner"];
			
			std::string currentSolverType = parameters_ ["current_solver_type"];
			std::string currentSolverMethod = parameters_ ["current_solver_method"];
			std::string currentSolverPreconditioner = parameters_ ["current_solver_preconditioner"];
			
			bool needSolverUpdate = 
				(desiredSolverType != currentSolverType) 
				|| 
				(desiredSolverMethod != currentSolverMethod)  
				||
				(desiredSolverPreconditioner != currentSolverPreconditioner);
			
			if (needSolverUpdate)
			{
				dolfin::log (dolfin::DBG, "Updating solver...");
				solver_ = createSolver ();
				dolfin::log (dolfin::DBG, "done!");
			}
			
			// define auxiliary string variables
			bool systemIsAssembled = parameters_ ["system_is_assembled"];
			bool forceReassembleSystem = parameters_ ["force_reassemble_system"];
			bool needReassemble = !systemIsAssembled || forceReassembleSystem;
			
			if (needReassemble)
			{
				dolfin::log (dolfin::DBG, "Reassembling system...");
				dolfin::assemble (*problemMatrix_, bilinearForm_);
				dolfin::assemble (rhsVector_, linearForm_);
				for (auto i : dirichletBCs_)
				{
					i.apply (*problemMatrix_, rhsVector_);
				}
				
				solver_ -> set_operator (problemMatrix_);
				
				parameters_ ["system_is_assembled"] = true;
				
				dolfin::log (dolfin::DBG, "done!");
			}
			
			dolfin::log (dolfin::DBG, "Solving problem...");
			solver_ -> solve (*solution_.vector (), rhsVector_);
			dolfin::log (dolfin::DBG, "done!");
		}



	template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
		void LinearDifferentialProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
		solve (const bool& mustReassemble)
		{
			parameters_ ["system_is_assembled"] = mustReassemble;
			solve ();
		}
	
	

	template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
		std::unique_ptr<dolfin::GenericLinearSolver> 
		LinearDifferentialProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
		createSolver ()
		{
			std::string solverType = parameters_["desired_solver_type"];
			std::string solverMethod = parameters_["desired_solver_method"];
			std::string solverPreconditioner = parameters_["desired_solver_preconditioner"];
			
			if (solverType == "lu_solver")
			{
				dolfin::log (dolfin::DBG, "Creating lu_solver...");
				std::unique_ptr<dolfin::GenericLinearSolver> solver (new dolfin::LUSolver (solverMethod));
				
				dolfin::log (dolfin::DBG, "Updating parameters...");
				parameters_ ["current_solver_type"] = solverType;
				parameters_ ["current_solver_method"] = solverMethod;
				parameters_ ["current_solver_preconditioner"] = solverPreconditioner;

				dolfin::log (dolfin::DBG, "done!");
				
				return solver;
			}
			else if (solverType == "krylov_solver")
			{
				dolfin::log (dolfin::DBG, "Creating krylov_solver...");
				std::unique_ptr<dolfin::GenericLinearSolver> solver (new dolfin::KrylovSolver (solverMethod, 
																							   solverPreconditioner));
				
				dolfin::log (dolfin::DBG, "Updating parameters...");
				parameters_ ["current_solver_type"] = solverType;
				parameters_ ["current_solver_method"] = solverMethod;
				parameters_ ["current_solver_preconditioner"] = solverPreconditioner;

				dolfin::log (dolfin::DBG, "done!");
				
				return solver;
			}
			else
			{
				dolfin::log (dolfin::DBG, "Creating solver of type %s...", solverType.c_str ());
				control_problem::LinearSolverFactory& factory = control_problem::LinearSolverFactory::Instance ();
				auto solver = factory.create (solverType);
				
				dolfin::log (dolfin::DBG, "Updating parameters...");
				parameters_ ["current_solver_type"] = solverType;
				parameters_ ["current_solver_method"] = solverMethod;
				parameters_ ["current_solver_preconditioner"] = solverPreconditioner;

				dolfin::log (dolfin::DBG, "done!");
				
				return solver;
			}
		}
}
#endif
