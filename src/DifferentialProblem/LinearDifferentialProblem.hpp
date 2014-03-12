#ifndef HH__LINEARDIFFERENTIALPROBLEM__HH
#define HH__LINEARDIFFERENTIALPROBLEM__HH

#include <dolfin.h>
#include <vector>
#include <string>
#include <memory>
#include <DifferentialProblem/AbstractDifferentialProblem.hpp>
#include <Factory/LinearSolverFactory.hpp>

namespace control_problem
{
	/*! \class LinearDifferentialProblem LinearDifferentialProblem.hpp
	 *  \brief Class for linear differential problems.
	 *
	 *  It inherits publicly from \c AbstractDifferentialProblem
	 *  and it extends its functionalities to a concrete differential
	 *  problem.
	 *  Template arguments are:
	 *  \arg T_BilinearForm the bilinear form type
	 *  \arg T_LinearForm the linear form type
	 *  This is needed because the dolfin function "coefficient_number" was neither declared
	 *  virtual nor implemented by FEniCS's developers in the dolfin class "Form", from which
	 *  we first thought to derive for the protected members of the LinearDifferentialProblem
	 *  class
	 */

	template <class T_BilinearForm, class T_LinearForm>
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
			 	 *  The stored mesh's and function space's ownership will be shared between the object and the input argument.
			 	 *  The bilinear and linear form will be created too, calling the constructor which takes the function space
			 	 *  as input.
			 	 */
				LinearDifferentialProblem (const std::shared_ptr<dolfin::Mesh> mesh, 
			                           	   const std::shared_ptr<dolfin::FunctionSpace> functionSpace);

				//! Constructor with references
				/*!
			 	 *  \param mesh the problem mesh as a const dolfin::Mesh&
			 	 *  \param functionSpace the problem finite element space as a const dolfin::FunctionSpace&
			 	 *  The stored mesh's and function space's ownership will be unique to the object, since the pointers are 
			 	 *  initialized using the \c new operator and mesh's and functionSpace's copy constructor
			 	 *  The bilinear and linear form will be created too, calling the constructor which takes the function space
			 	 *  as input.
			 	 */
				LinearDifferentialProblem (const dolfin::Mesh& mesh, 
			                           	   const dolfin::FunctionSpace& functionSpace);

				//! Constructor with rvalue references
				/*!
			 	 *  \param mesh the problem mesh as a dolfin::Mesh&&
			 	 *  \param functionSpace the problem finite element space as a dolfin::FunctionSpace&&
			 	 *  The stored mesh's and function space's ownership will be unique to the object, since the pointers are 
			 	 *  initialized using the \c new operator and mesh's and functionSpace's move constructor
			 	 *  The bilinear and linear form will be created too, calling the constructor which takes the function space
			 	 *  as input.
			 	 */
				LinearDifferentialProblem (dolfin::Mesh&& mesh, 
			                           	   dolfin::FunctionSpace&& functionSpace);

				//! Copy constructor
				/*! Note that this will copy the members that are shared pointers with the shared_ptr
			 	 * copy constructor. This means that, for example, the \c solver_ member will share
			 	 * ownership of the stored object with the contructor input argument.
				 * Use \c clone() method if you want a deep copy of the object:
				 * \code
				 *  LinearDifferentialProblem<T1, T2> foo; 
				 *  LinearDifferentialProblem<T1, T2> bar = foo.clone ();
				 * \endcode
			 	 */
				LinearDifferentialProblem (const LinearDifferentialProblem& rhs);

				//! Move constructor
				LinearDifferentialProblem (LinearDifferentialProblem&& rhs);


				/******************* GETTERS *******************/
				//! Get mesh as a const reference
				/*! 
			 	 *  \return a const reference to the problem's mesh
			 	 */
				const dolfin::Mesh& mesh () const;

				//! Get mesh as a shared_ptr
				/*! 
			 	 *  \return a const shared pointer to the problem's mesh
			 	 */
				const std::shared_ptr<const dolfin::Mesh> mesh () const;

				//! Get finite element space as a const reference
				/*! 
			 	 *  \return a const reference to the problem's function space
			 	 */
				const dolfin::FunctionSpace& functionSpace () const;

				//! Get finite element space as a shared_ptr
				/*! 
			 	 *  \return a const shared pointer to the problem's function space
			 	 */
				const std::shared_ptr<const dolfin::FunctionSpace> functionSpace () const;

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

				//! Set bilinear form
				/*! 
			 	 *  \param bilinearForm a const reference to a bilinear form on the
			 	 *                       problem's fnite element space 
			 	 */
				void setBilinearForm (const T_BilinearForm& bilinearForm);

				//! Set bilinear form [2]
				/*! 
			 	 *  \param bilinearForm a rvalue reference to a bilinear form on the
			 	 *                       problem's fnite element space 
			 	 */
				void setBilinearForm (T_BilinearForm&& bilinearForm);

				//! Set bilinear form coefficient [1]
				/*!
			 	 *  Set coefficient with given number. It simply wraps the set_coefficient dolfin function
			 	 *  \param i number of the coefficient to be set
			 	 *  \param coefficient coefficient's value
			 	 */
				void setBilinearFormCoefficient (std::size_t i, std::shared_ptr<const GenericFunction> coefficient);

				//! Set bilinear form coefficient [2]
				/*!
			 	 *  Set coefficient with given name. It simply wraps the set_coefficient dolfin function
			 	 *  \param name name of the coefficient to be set
			 	 *  \param coefficient coefficient's value
			 	 */
				void setBilinearFormCoefficient (std::string name, std::shared_ptr<const GenericFunction> coefficient);

				//! Set bilinear form coefficients
				/*!
			 	 *  Set coefficients in given map. It simply wraps the set_coefficient dolfin function
			 	 *  \param coefficients map of the coefficients to be set
			 	 */
				void setBilinearFormCoefficients (std::map<std::string, std::shared_ptr<const GenericFunction>> coefficients);

				//! Set linear form [1]
				/*! 
			 	 *  \param linearForm a const reference to a linear form on the
			 	 *                       problem's fnite element space 
			 	 */
				void setLinearForm (const T_LinearForm& linearForm);

				//! Set linear form [2]
				/*! 
			 	 *  \param linearForm a rvalue reference to a linear form on the
			 	 *                       problem's fnite element space 
			 	 */
				void setLinearForm (T_LinearForm&& linearForm);

				//! Set linear form coefficient [1]
				/*!
			 	 *  Set coefficient with given number. It simply wraps the set_coefficient dolfin function
			 	 *  \param i number of the coefficient to be set
			 	 *  \param coefficient coefficient's value
			 	 */
				void setLinearFormCoefficient (std::size_t i, std::shared_ptr<const GenericFunction> coefficient);

				//! Set linear form coefficient [2]
				/*!
			 	 *  Set coefficient with given name. It simply wraps the set_coefficient dolfin function
			 	 *  \param name name of the coefficient to be set
			 	 *  \param coefficient coefficient's value
			 	 */
				void setLinearFormCoefficient (std::string name, std::shared_ptr<const GenericFunction> coefficient);

				//! Set linear form coefficients
				/*!
			 	 *  Set coefficients in given map. It simply wraps the set_coefficient dolfin function
			 	 *  \param coefficients map of the coefficients to be set
			 	 */
				void setLinearFormCoefficients (std::map<std::string, std::shared_ptr<const GenericFunction>> coefficients);

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
				void removeDirichletBC (const std::vector<DirichletBC>::iterator& i);

				//! Set linear solver for the problem [1]
				/*!
			 	 *  \param solver a const reference to the solver to be set
			 	 */
				void setSolver (const dolfin::GenericLinearSolver& solver, std::string solverType);

				//! Set linear solver for the problem [2]
				/*!
			 	 *  \param solver a rvalue reference to the linear solver to be set
			 	 */
				void setSolver (dolfin::GenericLinearSolver&& solver, std::string solverType);

				//! Set linear solver for the problem [3]
				/*!
			 	 *  \param solver a rvalue reference to the linear solver to be set
			 	 */
				void setSolver (std::string solverType);

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
				
				//! Clone the current objecy
				/*!
				 *  This method will make a copy of the current object totally independent from it
				 *  \return a object of the same type as the one being copied
				 */
				LinearDifferentialProblem clone () const;

				//! Solve problem
				/*!
			 	 * This method solves the problem defined. It uses the private members' value to set the problem and then
			 	 * stores the solution in the private member \c solution_
			 	 */
				void solve ();

				//! Solve problem specifying flags
				/*!
			 	 *  \param mustReassemble true if the system operators (matrix and right hand side vector)
			 	 *         should be reassembled. It is false by default.
			 	 */
				void solve (const bool& mustReassemble);


				/******************* DESTRUCTOR *******************/

				virtual ~LinearDifferentialProblem () = default;

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

	template <class T_BilinearForm, class T_LinearForm>
		LinearDifferentialProblem<T_BilinearForm, T_LinearForm>::
		LinearDifferentialProblem (const std::shared_ptr<dolfin::Mesh> mesh, 
		                           const std::shared_ptr<dolfin::FunctionSpace> functionSpace) :
			mesh_ (mesh),
			functionSpace_ (functionSpace),
			bilinearForm_ (*functionSpace, *functionSpace),
			linearForm_ (*functionSpace),
			dirichletBCs_ (),
			solver_ (nullptr),
			solverType_ (),
			solution_ (functionSpace),
			isAssembled_ (0)
			problemMatrix_ (new dolfin::Matrix),
			rhsVector_ ()
		{ }



	template <class T_BilinearForm, class T_LinearForm>
		LinearDifferentialProblem<T_BilinearForm, T_LinearForm>::
		LinearDifferentialProblem (const dolfin::Mesh& mesh, 
		                           const dolfin::FunctionSpace& functionSpace) :
			mesh_ (new dolfin::Mesh (mesh)),
			functionSpace_ (new dolfin::FunctionSpace (functionSpace)),
			bilinearForm_ (*functionSpace, *functionSpace),
			linearForm_ (*functionSpace),
			dirichletBCs_ (),
			solver_ (nullptr),
			solverType_ (),
			solution_ (functionSpace),
			isAssembled_ (0)
			problemMatrix_ (new dolfin::Matrix),
			rhsVector_ ()
		{ }



	template <class T_BilinearForm, class T_LinearForm>
		LinearDifferentialProblem<T_BilinearForm, T_LinearForm>::
		LinearDifferentialProblem (dolfin::Mesh&& mesh, 
		                           dolfin::FunctionSpace&& functionSpace) :
			mesh_ (std::make_shared<dolfin::Mesh> (mesh)),
			functionSpace_ (std::make_shared<dolfin::FunctionSpace> (functionSpace)),
			bilinearForm_ (*functionSpace, *functionSpace),
			linearForm_ (*functionSpace),
			dirichletBCs_ (),
			solver_ (nullptr),
			solverType_ (),
			solution_ (functionSpace),
			isAssembled_ (0)
			problemMatrix_ (new dolfin::Matrix),
			rhsVector_ ()
		{ }



	template <class T_BilinearForm, class T_LinearForm>
		LinearDifferentialProblem<T_BilinearForm, T_LinearForm>::
		LinearDifferentialProblem (const LinearDifferentialProblem& rhs) : 
			mesh_ (rhs.mesh_),
			functionSpace_ (rhs.functionSpace_),
			bilinearForm_ (rhs.bilinearForm_),
			linearForm_ (rhs.linearForm_),
			dirichletBCs_ (rhs.dirichletBCs_),
			solver_ (nullptr),
			solverType_ (rhs.solverType_),
			solution_ (rhs.solution_),
			isAssembled_ (rhs.isAssembled_),
			problemMatrix_ (new dolfin::Matrix (rhs.problemMatrix_)),
			rhsVector_ (rhs.rhsVector_)
		{
			control_problem::LinearSolverFactory& factory = control_problem::LinearSolverFactory::Instance();
			solver_ = factory.create (rhs.solverType_);
		}



	template <class T_BilinearForm, class T_LinearForm>
		LinearDifferentialProblem<T_BilinearForm, T_LinearForm>::
		LinearDifferentialProblem (LinearDifferentialProblem&& rhs) :
			mesh_ (std::move (rhs.mesh_)),
			functionSpace_ (std::move (rhs.functionSpace_)),
			bilinearForm_ (std::move (rhs.bilinearForm_)),
			linearForm_ (std::move (rhs.linearForm_)),
			solver_ (std::move (rhs.solver_)),
			solverType_ (std::move (rhs.solverType_)),
			solution_ (std::move (rhs->solution_)),
			isAssembled_ (std::move (rhs.isAssembled_)),
			problemMatrix_ (std::move (rhs.problemMatrix_)),
			rhsVector_ (std::move (rhs.rhsVector_))
		{ }



	/******************* GETTERS *******************/

	template <class T_BilinearForm, class T_LinearForm>
		const dolfin::Mesh& LinearDifferentialProblem<T_BilinearForm, T_LinearForm>::mesh () const
		{
			return *mesh_;		
		}



	template <class T_BilinearForm, class T_LinearForm>
		const std::shared_ptr<const dolfin::Mesh> LinearDifferentialProblem<T_BilinearForm, T_LinearForm>::mesh () const
		{
			return mesh_;
		}



	template <class T_BilinearForm, class T_LinearForm>
		const dolfin::FunctionSpace& LinearDifferentialProblem<T_BilinearForm, T_LinearForm>::functionSpace () const
		{
			return *functionSpace_;
		}



	template <class T_BilinearForm, class T_LinearForm>
		const std::shared_ptr<const dolfin::FunctionSpace> LinearDifferentialProblem<T_BilinearForm, T_LinearForm>::
		functionSpace () const
		{
			return functionSpace_;
		}



	template <class T_BilinearForm, class T_LinearForm>
		const T_BilinearForm& LinearDifferentialProblem<T_BilinearForm, T_LinearForm>::bilinearForm () const
		{
			return bilinearForm_;
		}



	template <class T_BilinearForm, class T_LinearForm>
		const T_LinearForm& LinearDifferentialProblem<T_BilinearForm, T_LinearForm>::linearForm () const
		{
			return linearForm_;
		}



	template <class T_BilinearForm, class T_LinearForm>
		const dolfin::DirichletBC& LinearDifferentialProblem<T_BilinearForm, T_LinearForm>::
		dirichletBC (const std::size_t& i) const
		{
			return dirichletBCs_[i];
		}



	template <class T_BilinearForm, class T_LinearForm>
		const std::vector<const dolfin::DirichletBC>& LinearDifferentialProblem<T_BilinearForm, T_LinearForm>::
		dirichletBC () const
		{
			return dirichletBCs_;
		}



	template <class T_BilinearForm, class T_LinearForm>
		const dolfin::GenericLinearSolver& LinearDifferentialProblem<T_BilinearForm, T_LinearForm>::solver () const
		{
			return *solver_;
		}



	template <class T_BilinearForm, class T_LinearForm>
		const std::string& LinearDifferentialProblem<T_BilinearForm, T_LinearForm>::solverType () const
		{
			return solverType_;
		}

	template <class T_BilinearForm, class T_LinearForm>
		const dolfin::Function& LinearDifferentialProblem<T_BilinearForm, T_LinearForm>::solution () const
		{
			return solution_;
		}



	template <class T_BilinearForm, class T_LinearForm>
		const dolfin::Matrix& LinearDifferentialProblem<T_BilinearForm, T_LinearForm>::linearOperator () const
		{
			return *problemMatrix_;
		}



	template <class T_BilinearForm, class T_LinearForm>
		const dolfin::Vector& LinearDifferentialProblem<T_BilinearForm, T_LinearForm>::rhs () const
		{
			return rhsVector_;
		}



	/******************* SETTERS *******************/

	template <class T_BilinearForm, class T_LinearForm>
		bool LinearDifferentialProblem<T_BilinearForm, T_LinearForm>::setBilinearForm (const T_BilinearForm& bilinearForm)
		{
			bilinearForm_ = bilinearForm;
		}



	template <class T_BilinearForm, class T_LinearForm>
		bool LinearDifferentialProblem<T_BilinearForm, T_LinearForm>::setBilinearForm (T_BilinearForm&& bilinearForm)
		{
			bilinearForm_ = bilinearForm;
		}



	template <class T_BilinearForm, class T_LinearForm>
		void LinearDifferentialProblem<T_BilinearForm, T_LinearForm>::
		setBilinearFormCoefficient (std::size_t i, std::shared_ptr<const GenericFunction> coefficient)
		{
			bilinearForm_.set_coefficient (i, coefficient);
		}



	template <class T_BilinearForm, class T_LinearForm>
		void LinearDifferentialProblem<T_BilinearForm, T_LinearForm>::
		setBilinearFormCoefficient (std::string name, std::shared_ptr<const GenericFunction> coefficient)
		{
			bilinearForm_.set_coefficient (name, coefficient);
		}



	template <class T_BilinearForm, class T_LinearForm>
		void LinearDifferentialProblem<T_BilinearForm, T_LinearForm>::
		setBilinearFormCoefficients (std::map<std::string, std::shared_ptr<const GenericFunction>> coefficients)
		{
			bilinearForm_.set_coefficients (coefficients);
		}



	template <class T_BilinearForm, class T_LinearForm>
		void LinearDifferentialProblem<T_BilinearForm, T_LinearForm>::setLinearForm (const T_LinearForm& linearForm)
		{
			linearForm_ = linearForm;
		}



	template <class T_BilinearForm, class T_LinearForm>
		void LinearDifferentialProblem<T_BilinearForm, T_LinearForm>::setLinearForm (T_LinearForm&& linearForm)
		{
			linearForm_ = linearForm;
		}



	template <class T_BilinearForm, class T_LinearForm>
		void LinearDifferentialProblem<T_BilinearForm, T_LinearForm>::
		setLinearFormCoefficient (std::size_t i, std::shared_ptr<const GenericFunction> coefficient)
		{
			linearForm_.set_coefficient (i, coefficient);
		}



	template <class T_BilinearForm, class T_LinearForm>
		void LinearDifferentialProblem<T_BilinearForm, T_LinearForm>::
		setLinearFormCoefficient (std::string name, std::shared_ptr<const GenericFunction> coefficient)
		{
			linearForm_.set_coefficient (nome, coefficient);
		}



	template <class T_BilinearForm, class T_LinearForm>
		void LinearDifferentialProblem<T_BilinearForm, T_LinearForm>::
		setLinearFormCoefficients (std::map<std::string, std::shared_ptr<const GenericFunction>> coefficients)
		{
			linearForm_.set_coefficients (coefficients);
		}



	template <class T_BilinearForm, class T_LinearForm>
		void LinearDifferentialProblem<T_BilinearForm, T_LinearForm>::addDirichletBC (const dolfin::DirichletBC& dirichletCondition)
		{
			dirichletBCs_.emplace_back (dirichletCondition);
		}



	template <class T_BilinearForm, class T_LinearForm>
		void LinearDifferentialProblem<T_BilinearForm, T_LinearForm>::addDirichletBC (dolfin::DirichletBC&& dirichletCondition)
		{
			dirichletBCs_.emplace_back (dirichletCondition);
		}



	template <class T_BilinearForm, class T_LinearForm>
		void LinearDifferentialProblem<T_BilinearForm, T_LinearForm>::removeDirichletBC (const std::vector<DirichletBC>::iterator& i)
		{
			dirichletBCs_.erase (i);
		}



	template <class T_BilinearForm, class T_LinearForm>
		void LinearDifferentialProblem<T_BilinearForm, T_LinearForm>::
		setSolver (const dolfin::GenericLinearSolver& solver, std::string solverType)
		{
			solver_.reset (&solver);
			solverType_ = solverType;
		}



	template <class T_BilinearForm, class T_LinearForm>
		void LinearDifferentialProblem<T_BilinearForm, T_LinearForm>::
		setSolver (dolfin::GenericLinearSolver&& solver, std::solverType)
		{
			solver_.reset (&solver); 
			solverType_ = solverType;
		}


	
	template <class T_BilinearForm, class T_LinearForm>
		void setSolver (std::string solverType)
		{
			control_problem::LinearSolverFactory& factory = control_problem::LinearSolverFactory::Instance;
			solver_ = factory.create (rhs.solverType_);
		}
	


	template <class T_BilinearForm, class T_LinearForm>
		void LinearDifferentialProblem<T_BilinearForm, T_LinearForm>::
		setSolverParameters (const std::string& parameterName, const std::string& parameterValue)
		{
			solver_.parameters [parameterName] = parameterValue;	
		}



	template <class T_BilinearForm, class T_LinearForm>
		void LinearDifferentialProblem<T_BilinearForm, T_LinearForm>::
		setSolverParameters (const std::string& parameterName, const double& parameterValue)
		{
			solver_.parameters [parameterName] = parameterValue;	
		}



	template <class T_BilinearForm, class T_LinearForm>
		LinearDifferentialProblem<T_BilinearForm, T_LinearForm> LinearDifferentialProblem<T_BilinearForm, T_LinearForm>::
		clone () const
		{
			LinearDifferentialProblem newProblem (*mesh_, functionSpace_);
			
			newProblem.bilinearForm_ = bilinearForm_;
			newProblem.linearForm_   = linearForm_;
			newProblem.dirichletBCs_ = dirichletBCs_;
			
			control_problem::LinearSolverFactory& factory = control_problem::LinearSolverFactory::Instance();
			newProblem.solver_ = factory.create (solverType_);
			newProblem.solverType_ = solverType_;
			
			newProblem.solution_ = solution_;
			
			newProblem.isAssembled_ = isAssembled_;
		
			newProblem.problemMatrix_ = problemMatrix_;
			newProblem.rhsVector_ = rhsVector_;

			return newProblem;
		}



	template <class T_BilinearForm, class T_LinearForm>
		void LinearDifferentialProblem<T_BilinearForm, T_LinearForm>::solve () 
		{
			if (isAssembled_ == false)
			{
				dolfin::assemble (*problemMatrix_, bilinearForm_);
				dolfin::assemble (rhsVector_, linearForm_);
				for (auto i : dirichletBCs_)
				{
					i.apply (*problemMatrix_, rhsVector_);
				}
				isAssembled_ = true;
			}
			solver_ -> set_operator (problemMatrix_);
			solver_ -> (*solution_.vector (), rhsVector_);
		}



	template <class T_BilinearForm, class T_LinearForm>
		void LinearDifferentialProblem::<T_BilinearForm, T_LinearForm>::solve (const bool& mustReassemble)
		{
			isAssembled_ = false;
			this -> solve ();
		}
}
#endif
