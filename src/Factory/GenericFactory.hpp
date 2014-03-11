#ifndef HH__GENERICFACTORY__HH
#define HH__GENERICFACTORY__HH

#include <memory>
#include <functional>

namespace control_problem
{
	/*! class GenericFactory
	 * 
	 *  A generic template-ized factory class. 
	 *  Template arguments are
	 *  \arg the type of the object to be created
	 *  \arg the type of the identifier, that is the object that allows the factory to
	 *  distinguish among different possible objects to be built. For example, it could be a string containing
	 *  the name of the object
	 *  \arg the type of the builder. By default, this is set to a function that takes no arguments
	 *  and returns a unique_ptr to the created object
	 * 
	 * ciao
	 */
	template 
		<
			T_AbstractProduct, 
			T_Identifier, 
			T_Builder = std::function <std::unique_ptr<T_AbstractProduct> ()> 
		>
		class GenericFactory
		{
			private:
				GenericFactory () = default;
				GenericFactory (GenericFactory&& factory) = delete;
				GenericFactory ()<++>
		}<++>
}
#endif
