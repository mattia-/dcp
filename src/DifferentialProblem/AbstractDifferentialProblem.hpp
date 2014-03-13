#ifndef HH__ABSTRACTDIFFERENTIALPROBLEM__HH
#define HH__ABSTRACTDIFFERENTIALPROBLEM__HH

#include <memory>
#include <dolfin.h>

namespace control_problem
{
	/*! \class AbstractDifferentialProblem AbstractDifferentialProblem.hpp
	 *  \brief Abstract base class for differential problems. 
	 *         
	 *  This class contains the basic interface for a differential problem to be
	 *  solved with FEniCS library. It is an abstract class, it only provides the
	 *  basic interface to all differential problems
	 */ 
	class AbstractDifferentialProblem
	{
		// ---------------------------------------------------------------------------------------------//

		public:
			//! Default constructor
			AbstractDifferentialProblem () = default;
			
			//! Constructor that takes a function space as input argument
			AbstractDifferentialProblem (const dolfin::FunctionSpace& functionSpace);
			
			//! Constrtuctor that takes a function as input argument [1] (lvalue reference version)
			AbstractDifferentialProblem (const dolfin::Function& solution);
			
			//! Constrtuctor that takes a function as input argument [2] (rvalue reference version)
			AbstractDifferentialProblem (dolfin::Function&& solution);

			//! Destructor
			/*! Default destructor, since members of the class are trivially 
			 * destructible.
			 * It is declared virtual so that derived classes' constructor
			 * can be called on derived classes.
			 * The "default-ness" is set in implementation outside of the class for compatibility with
			 * gcc-4.6, which does not allow virtual members to be defaulted in class
			 */
			virtual ~AbstractDifferentialProblem ();

			//! Solve method
			/*!
			 * Solves differential problem storing the solution 
			 * in the private member solution.
			 * It is a pure virtual method that needs to be overridden
			 * in any concrete instance of the class
			 */
			virtual void solve () = 0;

			// ---------------------------------------------------------------------------------------------//

		protected:
			//! Solution of the differential problem
			dolfin::Function solution_;

			// ---------------------------------------------------------------------------------------------//

		private:

	};
	
	// =============================================================== //
	// ====================== IMPLEMENTATION ========================= //
	// =============================================================== //
	
	// ***** CONSTRUCTORS
	AbstractDifferentialProblem::AbstractDifferentialProblem (const dolfin::FunctionSpace& functionSpace) : 
		solution_ (functionSpace)
	{  }



	AbstractDifferentialProblem::AbstractDifferentialProblem (const dolfin::Function& solution) :
		solution_ (solution)
	{  }



	AbstractDifferentialProblem::AbstractDifferentialProblem (dolfin::Function&& solution) :
		solution_ (std::move (solution))
	{  }
	
	
	
	// ***** DESTRUCTORS
	// this is done for compatibility with gcc 4.6, which doesn't allow virtual members to be defualted in class body
	AbstractDifferentialProblem::~AbstractDifferentialProblem () = default;
}
#endif
