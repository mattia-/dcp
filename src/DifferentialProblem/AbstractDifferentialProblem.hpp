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
			/*! 
			 * Default constructor provided by \c C++
			 */
			AbstractDifferentialProblem () = default;

			//! Destructor
			/*! Empty destructor, since members of the class are trivially 
			 * destructible.
			 * It is declared virtual so that derived classes' constructor
			 * can be called on derived classes
			 */
			virtual ~AbstractDifferentialProblem () = default;

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
}
#endif
