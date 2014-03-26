#ifndef SRC_UTILS_COMPARISONFUNCTORS_HPP_INCLUDE_GUARD
#define SRC_UTILS_COMPARISONFUNCTORS_HPP_INCLUDE_GUARD

/*! This file contains definitions of the function object \c control_problem::less, 
 * \c control_problem::greater, \c control_problem::equal_to, ... so that we can use them in 
 * standard library containers that need an ordering relationship throughout this
 * project. We prefer not to specialize comparison operators in std library (e.g. std::less<T>,
 * std::greater<T>, ...) because otherwise projects that use this library will also inherit 
 * such specializations
 */
#include <functional>
#include <utility>
#include <string>

namespace control_problem
{
	/*! \class less
	 *  \brief class that implements the "less than" comparison. 
	 *  
	 *  Its call operator returns true if first argument is less than the second argument, 0 otherwise
	 *  This is not a template class. Instead, call operator is overloaded for every type we need in this
	 *  project. If you want more generality, use \c std::less<T>
	 */
	class less
	{
		public:
			//! Comparison operator for objects of type \c std::tuple<std::string, std::string, std::string>
			bool operator() (const std::tuple <std::string, std::string, std::string>& lhs, 
			                 const std::tuple <std::string, std::string, std::string>& rhs);
	};
}
#endif
