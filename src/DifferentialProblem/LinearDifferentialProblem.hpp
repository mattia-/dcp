#ifndef HH__LINEARDIFFERENTIALPROBLEM__HH
#define HH__LINEARDIFFERENTIALPROBLEM__HH

#include <dolfin.h>
#include <vector>
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

				//!  Constructor with shared pointers
				/*!
			 	 *  \param mesh the problem mesh as a const std::shared_ptr to dolfin::Mesh
			 	 *  \param functionSpace the problem finite element space as a const std::shared_ptr to dolfin::FunctionSpace
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
			                           	   const std::string& solverType,
			                           	   const std::string& solverMethod,
			                           	   const std::string& solverPreconditioner);
				

				//! Constructor with references
				/*!
			 	 *  \param mesh the problem mesh as a const dolfin::Mesh&
			 	 *  \param functionSpace the problem finite element space as a const dolfin::FunctionSpace&
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
			                           	   const std::string& solverType,
			                           	   const std::string& solverMethod,
			                           	   const std::string& solverPreconditioner);

				//! Constructor with rvalue references
				/*!
			 	 *  \param mesh the problem mesh as a dolfin::Mesh&&
			 	 *  \param functionSpace the problem finite element space as a dolfin::FunctionSpace&&
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
				const T_BilinearForm& bilinearForm () const;

				//! Get const reference to the problem's linear form
				/*! 
			 	 *  \return a const reference to the problem's linear form
			 	 */
				const T_LinearForm& linearForm () const;

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

				//! Get const reference to the problem's solver 
				/*! 
			 	 *  \return a const reference to the problem's solver.
			 	 */
				const dolfin::GenericLinearSolver& solver () const;

				//! Get const reference to the problem's solver's type
				/*! 
			 	 *  \return a const reference to the problem's solver's type.
			 	 */
				const std::string& solverType () const;

				//! Get const reference to the problem's solution
				/*!
			 	 *  \return a const reference to the problem's solution
			 	 */
				const dolfin::Function& solution () const; 	

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

				//! Set bilinear form coefficient [1]
				/*!
			 	 *  Set coefficient with given number. It simply wraps the set_coefficient dolfin function
			 	 *  \param i number of the coefficient to be set
			 	 *  \param coefficient coefficient's value. We use boost::shared_ptr instead of std::shared_ptr for
				 *  compatibility with dolfin1.3.0
			 	 */
				void setBilinearFormCoefficient (const std::size_t& i, 
				                                 const boost::shared_ptr<const dolfin::GenericFunction> coefficient);

				//! Set bilinear form coefficient [2]
				/*!
			 	 *  Set coefficient with given name. It simply wraps the set_coefficient dolfin function
			 	 *  \param name name of the coefficient to be set
			 	 *  \param coefficient coefficient's value. We use boost::shared_ptr instead of std::shared_ptr for
				 *  compatibility with dolfin1.3.0
			 	 */
				void setBilinearFormCoefficient (const std::string& name, 
				                                 const boost::shared_ptr<const dolfin::GenericFunction> coefficient);

				//! Set bilinear form coefficients
				/*!
			 	 *  Set coefficients in given map. It simply wraps the set_coefficient dolfin function
			 	 *  \param coefficients map of the coefficients to be set. 
				 *  We use boost::shared_ptr instead of std::shared_ptr for compatibility with dolfin1.3.0
			 	 */
				void setBilinearFormCoefficients 
					(const std::map<std::string, boost::shared_ptr<const dolfin::GenericFunction>>& coefficients);

				//! Set linear form coefficient [1]
				/*!
			 	 *  Set coefficient with given number. It simply wraps the set_coefficient dolfin function
			 	 *  \param i number of the coefficient to be set
			 	 *  \param coefficient coefficient's value. We use boost::shared_ptr instead of std::shared_ptr 
				 *  for compatibility with dolfin1.3.0
			 	 */
				void setLinearFormCoefficient (const std::size_t& i, 
				                               const boost::shared_ptr<const dolfin::GenericFunction> coefficient);

				//! Set linear form coefficient [2]
				/*!
			 	 *  Set coefficient with given name. It simply wraps the set_coefficient dolfin function
			 	 *  \param name name of the coefficient to be set
			 	 *  \param coefficient coefficient's value.
				 *  We use boost::shared_ptr instead of std::shared_ptr for compatibility with dolfin1.3.0
			 	 */
				void setLinearFormCoefficient (const std::string& name, 
				                               const boost::shared_ptr<const dolfin::GenericFunction> coefficient);

				//! Set linear form coefficients
				/*!
			 	 *  Set coefficients in given map. It simply wraps the set_coefficient dolfin function
			 	 *  \param coefficients map of the coefficients to be set.
				 *  We use boost::shared_ptr instead of std::shared_ptr for compatibility with dolfin1.3.0
			 	 */
				void setLinearFormCoefficients 
					(const std::map<std::string, boost::shared_ptr<const dolfin::GenericFunction>>& coefficients);
				
				//! Set integration subdomains for the bilinear form
				/*! Input arguments are:
				 *  \param meshFunction the mesh function used to set the integration subdomains. Template-ized over
				 *  the type of stored objects
				 *  \param subdomainType the type of the subdomains, chosen among those provided by the enumeration
				 *  class \c control_problem::SubdomainType
				 */
				template <class T>
					void setBilinearFormIntegrationSubdomains (const dolfin::MeshFunction<T>& meshFunction, 
					                                           const control_problem::SubdomainType& subdomainType);

				//! Set integration subdomains for the linear form
				/*! Input arguments are:
				 *  \param meshFunction the mesh function used to set the integration subdomains. Template-ized over
				 *  the type of stored objects
				 *  \param subdomainType the type of the subdomains, chosen among those provided by the enumeration
				 *  class \c SubdomainType
				 */
				template <class T>
					void setLinearFormIntegrationSubdomains (const dolfin::MeshFunction<T>& meshFunction, 
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

				//! Set linear solver for the problem
				/*!
			 	 *  \param solver a rvalue reference to the linear solver to be set
			 	 */
				void setSolver (const std::string& solverType, 
				                const std::string& solverMethod,
				                const std::string& solverPreconditioner);

				//! Set linear solver parameters [1]
				/*!
			 	 *  \param parameterName string containing the name of the linear solver parameter to be set
			 	 *  \param parameterValue string containing the value of the linear solver parameter to be set
			 	 */
				void setSolverParameters (const std::string& parameterName, const std::string& parameterValue);

				//! Set linear solver parameters [2]
				/*!
			 	 *  \param parameterName string containing the name of the linear solver parameter to be set
			 	 *  \param parameterValue double containing the value of the linear solver parameter to be set
			 	 */
				void setSolverParameters (const std::string& parameterName, const double& parameterValue);

				
				/******************* METHODS *******************/

				//! Solve problem
				/*!
			 	 * This method solves the problem defined. It uses the private members' value to set the problem and then
			 	 * stores the solution in the private member \c solution_
			 	 */
				void solve ();

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
				 *  \param solverType the type of the solver. If not given, the default is lu_solver
				 *  \param solverMethod the method to be used. The default value will create a default (in the sense 
				 *  of dolfin library) solver of the type passed as first argument
				 *  \param solverPreconditioner the preconditioner to be used.
				 *  Note that if solver type is neither lu_solver nor krylov_solver, the solver type will be looked for
				 *  in the linear solver factory of type passed to the class as template argument, and the \c create
				 *  method of such factory will be called. 
				 */
				std::unique_ptr<dolfin::GenericLinearSolver> createSolver (const std::string& solverType, 
				                                                           const std::string& solverMethod,
				                                                           const std::string& solverPreconditioner);
				
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
				T_BilinearForm bilinearForm_;

				//! The linear form
				T_LinearForm linearForm_;

				//! The Dirichlet's boundary conditions vector
				std::vector<dolfin::DirichletBC> dirichletBCs_;

				//! The solver
				/*! 
			 	 *  We use a pointer so that polymorphism can be applied.
			 	 */
				std::unique_ptr<dolfin::GenericLinearSolver> solver_;
				
				//! The solver's type, as a string
				std::string solverType_;

				//! Boolean flag to indicate if system has already been assembled
				/*! 
			 	 *  It is equal to 1 if the system has already been assembled, 0 otherwise.
			 	 *  Note that the system will not be assembled until solve () is called for
			 	 *  the first time.
			 	 */
				bool isAssembled_;

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



	// ==============================================================================================//
	// ==================================== IMPLEMENTATION ==========================================//
	// ==============================================================================================//


	/******************* CONSTRUCTORS *******************/

	template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
		LinearDifferentialProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
		LinearDifferentialProblem (const std::shared_ptr<dolfin::Mesh> mesh, 
		                           const std::shared_ptr<dolfin::FunctionSpace> functionSpace,
			                       const std::string& solverType = "lu_solver",
			                       const std::string& solverMethod = "default",
			                       const std::string& solverPreconditioner = "default") :
			AbstractDifferentialProblem (functionSpace),
			mesh_ (mesh),
			functionSpace_ (functionSpace),
			bilinearForm_ (*functionSpace, *functionSpace),
			linearForm_ (*functionSpace),
			dirichletBCs_ (),
			solver_ (createSolver (solverType, solverMethod, solverPreconditioner)),
			solverType_ (solverType),
			isAssembled_ (0),
			problemMatrix_ (new dolfin::Matrix),
			rhsVector_ ()
		{ }



	template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
		LinearDifferentialProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
		LinearDifferentialProblem (const dolfin::Mesh& mesh, 
		                           const dolfin::FunctionSpace& functionSpace,
			                       const std::string& solverType = "lu_solver",
			                       const std::string& solverMethod = "default",
			                       const std::string& solverPreconditioner = "default") :
			AbstractDifferentialProblem (functionSpace),
			mesh_ (new dolfin::Mesh (mesh)),
			functionSpace_ (new dolfin::FunctionSpace (functionSpace)),
			bilinearForm_ (functionSpace, functionSpace),
			linearForm_ (functionSpace),
			dirichletBCs_ (),
			solver_ (createSolver (solverType, solverMethod, solverPreconditioner)),
			solverType_ (solverType),
			isAssembled_ (0),
			problemMatrix_ (new dolfin::Matrix),
			rhsVector_ ()
		{ }



	template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
		LinearDifferentialProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
		LinearDifferentialProblem (dolfin::Mesh&& mesh, 
		                           dolfin::FunctionSpace&& functionSpace,
			                       const std::string& solverType = "lu_solver",
			                       const std::string& solverMethod = "default",
			                       const std::string& solverPreconditioner = "default") :
			AbstractDifferentialProblem (functionSpace),
			mesh_ (std::make_shared<dolfin::Mesh> (mesh)),
			functionSpace_ (std::make_shared<dolfin::FunctionSpace> (functionSpace)),
			bilinearForm_ (functionSpace, functionSpace),
			linearForm_ (functionSpace),
			dirichletBCs_ (),
			solver_ (createSolver (solverType, solverMethod, solverPreconditioner)),
			solverType_ (solverType),
			isAssembled_ (0),
			problemMatrix_ (new dolfin::Matrix),
			rhsVector_ ()
		{ }



	template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
		LinearDifferentialProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
		LinearDifferentialProblem (const std::shared_ptr<dolfin::Mesh> mesh, 
		                           const std::shared_ptr<dolfin::FunctionSpace> functionSpace,
		                           const T_BilinearForm& bilinearForm,
		                           const T_LinearForm& linearForm,
		                           const std::string& solverType = "lu_solver",
		                           const std::string& solverMethod = "default",
		                           const std::string& solverPreconditioner = "default") :
			AbstractDifferentialProblem (functionSpace),
			mesh_ (mesh),
			functionSpace_ (functionSpace),
			bilinearForm_ (bilinearForm),
			linearForm_ (linearForm),
			dirichletBCs_ (),
			solver_ (createSolver (solverType, solverMethod, solverPreconditioner)),
			solverType_ (solverType),
			isAssembled_ (0),
			problemMatrix_ (new dolfin::Matrix),
			rhsVector_ ()
	{ }



	template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
		LinearDifferentialProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
		LinearDifferentialProblem (const dolfin::Mesh& mesh, 
		                           const dolfin::FunctionSpace& functionSpace,
		                           const T_BilinearForm& bilinearForm,
		                           const T_LinearForm& linearForm,
			                       const std::string& solverType = "lu_solver",
			                       const std::string& solverMethod = "default",
			                       const std::string& solverPreconditioner = "default") :
			AbstractDifferentialProblem (functionSpace),
			mesh_ (new dolfin::Mesh (mesh)),
			functionSpace_ (new dolfin::FunctionSpace (functionSpace)),
			bilinearForm_ (bilinearForm),
			linearForm_ (linearForm),
			dirichletBCs_ (),
			solver_ (createSolver (solverType, solverMethod, solverPreconditioner)),
			solverType_ (solverType),
			isAssembled_ (0),
			problemMatrix_ (new dolfin::Matrix),
			rhsVector_ ()
	{ }



	template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
		LinearDifferentialProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
		LinearDifferentialProblem (dolfin::Mesh&& mesh, 
		                           dolfin::FunctionSpace&& functionSpace,
		                           T_BilinearForm&& bilinearForm,
		                           T_LinearForm&& linearForm,
		                           const std::string& solverType = "lu_solver",
		                           const std::string& solverMethod = "default",
		                           const std::string& solverPreconditioner = "default") :
			AbstractDifferentialProblem (functionSpace),
			mesh_ (std::make_shared<dolfin::Mesh> (mesh)),
			functionSpace_ (std::make_shared<dolfin::FunctionSpace> (functionSpace)),
			bilinearForm_ (std::move (bilinearForm)),
			linearForm_ (std::move (linearForm)),
			dirichletBCs_ (),
			isAssembled_ (0),
			solver_ (createSolver (solverType, solverMethod, solverPreconditioner)),
			solverType_ (solverType),
			problemMatrix_ (new dolfin::Matrix),
			rhsVector_ ()
	{ }

	
	/******************* DESTRUCTOR *******************/

	// this is done for compatibility with gcc-4.6, which doesn't allow virtual members to be defualted in class body
	template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
		LinearDifferentialProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
		~LinearDifferentialProblem () = default;
	
	


	/******************* GETTERS *******************/

	template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
		const dolfin::Mesh& LinearDifferentialProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
		mesh () const
		{
			return *mesh_;		
		}



	template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
		const dolfin::FunctionSpace& LinearDifferentialProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
		functionSpace () const
		{
			return *functionSpace_;
		}



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
		const dolfin::DirichletBC& LinearDifferentialProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
		dirichletBC (const std::size_t& i) const
		{
			return dirichletBCs_[i];
		}



	template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
		const std::vector<const dolfin::DirichletBC>& 
		LinearDifferentialProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
		dirichletBC () const
		{
			return dirichletBCs_;
		}



	template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
		const dolfin::GenericLinearSolver& LinearDifferentialProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
		solver () const
		{
			return *solver_;
		}



	template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
		const std::string& LinearDifferentialProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
		solverType () const
		{
			return solverType_;
		}

	template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
		const dolfin::Function& LinearDifferentialProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
		solution () const
		{
			return solution_;
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
		setBilinearFormCoefficient (const std::size_t& i, const boost::shared_ptr<const dolfin::GenericFunction> coefficient)
		{
			bilinearForm_.set_coefficient (i, coefficient);
		}



	template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
		void LinearDifferentialProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
		setBilinearFormCoefficient (const std::string& name, const boost::shared_ptr<const dolfin::GenericFunction> coefficient)
		{
			bilinearForm_.set_coefficient (name, coefficient);
		}



	template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
		void LinearDifferentialProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
		setBilinearFormCoefficients (const std::map<std::string, boost::shared_ptr<const dolfin::GenericFunction>>& coefficients)
		{
			bilinearForm_.set_coefficients (coefficients);
		}



	template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
		void LinearDifferentialProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
		setLinearFormCoefficient (const std::size_t& i, const boost::shared_ptr<const dolfin::GenericFunction> coefficient)
		{
			linearForm_.set_coefficient (i, coefficient);
		}



	template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
		void LinearDifferentialProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
		setLinearFormCoefficient (const std::string& name, const boost::shared_ptr<const dolfin::GenericFunction> coefficient)
		{
			linearForm_.set_coefficient (name, coefficient);
		}



	template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
		void LinearDifferentialProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
		setLinearFormCoefficients (const std::map<std::string, boost::shared_ptr<const dolfin::GenericFunction>>& coefficients)
		{
			linearForm_.set_coefficients (coefficients);
		}



	template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
		template <class T>
			void LinearDifferentialProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
			setBilinearFormIntegrationSubdomains (const dolfin::MeshFunction<T>& meshFunction, 
			                                      const control_problem::SubdomainType& subdomainType)
			{
				if (subdomainType == control_problem::SubdomainType::INTERNAL_CELLS)
				{
					bilinearForm_.dx = meshFunction;
				}
				else if (subdomainType == control_problem::SubdomainType::INTERNAL_FACETS)
				{
					bilinearForm_.dS = meshFunction;
				}
				else if (subdomainType == control_problem::SubdomainType::BOUNDARY_FACETS)
				{
					bilinearForm_.ds = meshFunction;
				}
				else
				{
					std::cerr 
						<< "Warning: unknown subdomain type when trying to apply mesh function to bilinear form" 
						<< std::endl;
				}
			}



	template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
		template <class T>
			void LinearDifferentialProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
			setLinearFormIntegrationSubdomains (const dolfin::MeshFunction<T>& meshFunction, 
		                                    	const control_problem::SubdomainType& subdomainType)
			{
				if (subdomainType == control_problem::SubdomainType::INTERNAL_CELLS)
					{
						linearForm_.dx = meshFunction;
					}
					else if (subdomainType == control_problem::SubdomainType::INTERNAL_FACETS)
					{
						linearForm_.dS = meshFunction;
					}
					else if (subdomainType == control_problem::SubdomainType::BOUNDARY_FACETS)
					{
						linearForm_.ds = meshFunction;
					}
					else
					{
						std::cerr 
							<< "Warning: unknown subdomain type when trying to apply mesh function to linear form" 
							<< std::endl;
					}
				}



	template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
		void LinearDifferentialProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
		addDirichletBC (const dolfin::DirichletBC& dirichletCondition)
		{
			dirichletBCs_.emplace_back (dirichletCondition);
		}



	template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
		void LinearDifferentialProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
		addDirichletBC (dolfin::DirichletBC&& dirichletCondition)
		{
			dirichletBCs_.emplace_back (dirichletCondition);
		}



	template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
		void LinearDifferentialProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
		removeDirichletBC (const std::vector<dolfin::DirichletBC>::iterator& i)
		{
			dirichletBCs_.erase (i);
		}



	template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
		void LinearDifferentialProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
		setSolver (const std::string& solverType = "lu_solver",
		           const std::string& solverMethod = "default",
		           const std::string& solverPreconditioner = "default")
		{
			solver_ = createSolver (solverType, solverMethod, solverPreconditioner);
			solverType_ = solverType;
		}
	


	template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
		void LinearDifferentialProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
		setSolverParameters (const std::string& parameterName, const std::string& parameterValue)
		{
			solver_ -> parameters [parameterName] = parameterValue;	
		}



	template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
		void LinearDifferentialProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
		setSolverParameters (const std::string& parameterName, const double& parameterValue)
		{
			solver_ -> parameters [parameterName] = parameterValue;	
		}



	template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
		void LinearDifferentialProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
		solve () 
		{
			if (isAssembled_ == false)
			{
				dolfin::assemble (*problemMatrix_, bilinearForm_);
				dolfin::assemble (rhsVector_, linearForm_);
				for (auto i : dirichletBCs_)
				{
					i.apply (*problemMatrix_, rhsVector_);
				}
				
				solver_ -> set_operator (problemMatrix_);
				
				isAssembled_ = true;
			}
			solver_ -> solve (*solution_.vector (), rhsVector_);
		}



	template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
		void LinearDifferentialProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
		solve (const bool& mustReassemble)
		{
			isAssembled_ = false;
			solve ();
		}
	

	template <class T_BilinearForm, class T_LinearForm, class T_LinearSolverFactory>
		std::unique_ptr<dolfin::GenericLinearSolver> 
		LinearDifferentialProblem<T_BilinearForm, T_LinearForm, T_LinearSolverFactory>::
		createSolver (const std::string& solverType = "lu_solver", 
		              const std::string& solverMethod = "default", 
		              const std::string& solverPreconditioner = "default")
		{
			if (solverType == "lu_solver")
			{
				return std::unique_ptr<dolfin::GenericLinearSolver> (new dolfin::LUSolver (solverMethod));
			}
			else if (solverType == "krylov_solver")
			{
				return std::unique_ptr<dolfin::GenericLinearSolver> 
					(new dolfin::KrylovSolver (solverMethod, solverPreconditioner));
			}
			else
			{
				control_problem::LinearSolverFactory& factory = control_problem::LinearSolverFactory::Instance ();
				return factory.create (solverType);
			}
		}
}
#endif
