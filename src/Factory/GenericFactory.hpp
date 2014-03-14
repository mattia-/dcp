#ifndef HH__GENERICFACTORY__HH
#define HH__GENERICFACTORY__HH

#include <memory>
#include <functional>
#include <map>
#include <vector>
#include <iostream>

namespace control_problem
{
	/*! class GenericFactory
	 *  \brief class containing a generic factory
	 * 
	 *  A generic template-ized factory class. 
	 *  Template arguments are
	 *  \arg the type of the object to be created
	 *  \arg the type of the identifier, that is the object that allows the factory to
	 *  distinguish among different possible objects to be built. For example, it could be a string containing
	 *  the name of the object
	 *  \arg the type of the builder. By default, this is set to a function that takes no arguments
	 *  and returns a unique_ptr to the created object
	 */
	template 
		<
			class T_AbstractProduct_, 
			class T_Identifier_, 
			class T_Builder_ = std::function <std::unique_ptr<T_AbstractProduct_> ()> 
		>
			class GenericFactory
			{
				public:
					typedef T_AbstractProduct_ T_AbstractProduct;
					typedef T_Identifier_      T_Identifier;
					typedef T_Builder_         T_Builder;

					//! Method to access the only instance of the factory
					static GenericFactory& Instance();

					//! Get the object with given identifier 
					/*!
					  The pointer is null if no match was found for the identifier.
					  */
					std::unique_ptr<T_AbstractProduct> create (const T_Identifier& identifier) const;

					//! Register the given rule
					void add (const T_Identifier& identifier, const T_Builder& builder);

					//! Returns a list of registered rules
					std::vector<T_Identifier> registered () const;

					//! Unregister a rule
					void unregister (const T_Identifier& name);

					//! Destructor
					~GenericFactory() = default;

				private:
					//! Default empty constructor
					/*! Declared private because we want only methods of the class to be able to use it.
					 *  In particular, we will use it in the Instance() method.
					 *  This is needed because we want the class to be a singleton
					 */
					GenericFactory () = default;

					//! Copy constructor is deleted, since the factory is a singleton
					GenericFactory (const GenericFactory& factory) = delete;

					//! Move constructor is deleted, since the factory is a singleton
					GenericFactory (GenericFactory&& factory) = delete;

					//! Copy operator deleted, since the factory it is a Singleton
					GenericFactory& operator= (const GenericFactory& factory) = delete;

					//! It contains the actual object factory.
					std::map<T_Identifier, T_Builder> storedData_;
			};

	// ======================================================================= //
	// ============================== IMPLEMENTATION ========================= //
	// ======================================================================= //


	// We use Meyer's trick to istantiate the factory.
	template
		<
			class T_AbstractProduct,
			class T_Identifier,
			class T_Builder
		>
		GenericFactory<T_AbstractProduct, T_Identifier, T_Builder>& 
		GenericFactory<T_AbstractProduct, T_Identifier, T_Builder>::Instance() 
		{
			static GenericFactory factory;
			return factory;
		}
	


	template
		<
			class T_AbstractProduct,
			class T_Identifier,
			class T_Builder
		>
		std::unique_ptr <T_AbstractProduct> 
		GenericFactory<T_AbstractProduct, T_Identifier, T_Builder>::create(const T_Identifier& identifier) const
		{
			auto f = storedData_.find (identifier); 
			if (f == storedData_.end ())
			{
				std::cerr << "Warning: Identifier " << identifier << " not found in factory" << std::endl;
				return nullptr;
			}
			else
			{
				return std::unique_ptr<T_AbstractProduct> (f->second());
			}
		}
	


	template
		<
			class T_AbstractProduct,
			class T_Identifier,
			class T_Builder
		>
		void 
		GenericFactory<T_AbstractProduct, T_Identifier, T_Builder>::
		add (const T_Identifier& identifier, const T_Builder& builder)
		{
			auto f = storedData_.insert(std::make_pair(identifier, builder));
		    if (f.second == false)
			{
				std::cerr 
					<< "Warning: double registration in factory while trying to register product "
					<< std::endl;
			}
		}


	
	template
		<
			typename T_AbstractProduct,
			typename T_Identifier,
			typename T_Builder
		>
		std::vector<T_Identifier> 
		GenericFactory<T_AbstractProduct, T_Identifier, T_Builder>::registered () const
		{
			std::vector<T_Identifier> tmp;
			tmp.reserve (storedData_.size());
			
			for(auto i = storedData_.begin(); i != storedData_.end(); ++i)
			    tmp.push_back (i->first);
			
			return tmp;
		}
};

#endif
