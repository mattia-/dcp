#ifndef HH__PROXY__HH
#define HH__PROXY__HH

#include "GenericFactory.hpp"

namespace control_problem
{
	/*! \class Proxy
	 *  \brief A class that provides registration into a factory.
	 *  
	 *  This class provides an easy and quick way to register 
	 *  objects in a generic factory, independent of the classes
	 *  used to instantiate the factory. 
	 *  To register an object into a factory, one only needs to 
	 *  create a proxy object feeding the constructor with the
	 *  builder one wants to use for the creation of such object.
	 */
	template <class T_Factory, class T_ConcreteProduct>
		class Proxy 
		{
			public:
				//! The constructor, which also does the registration
				Proxy (const T_Identifier_type& identifier);

				//! The builder. It must comply with the signature which we
				// want to use for our factory
				static std::unique_ptr<T_Factory::T_AbstractProduct> Build()
				{
					return std::unique_ptr<T_Factory::T_AbstractProducte> (new T_ConcreteProduct ());
				}

			private:
				Proxy (const Proxy& proxy) = delete;
				Proxy (Proxy&& proxy) = delete;
				Proxy& operator= (const Proxy& proxy) = delete;
		};


	template <typename T_Factory, typename T_ConcreteProduct>
		Proxy<T_Factory,T_ConcreteProduct>::Proxy (const T_Factory::T_Identifier_type& identifier) 
		{
			// get the factory. First time creates it.
			T_Factory& factory (T_Factory::Instance ());
			
			// Insert the builder.
			factory.add (identifier, Proxy <T_Factory,T_ConcreteProduct>::Build);
		}
};

#endif
