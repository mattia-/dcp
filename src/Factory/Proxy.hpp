#ifndef SRC_FACTORY_PROXY_HPP_INCLUDE_GUARD
#define SRC_FACTORY_PROXY_HPP_INCLUDE_GUARD

#include <Factory/GenericFactory.hpp>

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
     *  Template arguments are:
     *  \arg the type of the factory
     *  \arg the type of the concrete product
     */
    template <class T_Factory_, class T_ConcreteProduct_>
        class Proxy 
        {
            public:
                typedef typename T_Factory_::T_AbstractProduct T_Factory_AbstractProduct;
                typedef typename T_Factory_::T_Identifier      T_Factory_Identifier;
                typedef typename T_Factory_::T_Builder         T_Factory_Builder;
                typedef          T_Factory_                    T_Factory;
                typedef          T_ConcreteProduct_            T_ConcreteProduct;
                
                //! The constructor, which also does the registration
                /*! Input arguments are:
                 *  \param the identifier that will identify the object in the factory
                 */
                Proxy (const T_Factory_Identifier& identifier);

                //! The builder for the default construction. It must comply with the signature which we
                //! want to use for our factory
                static std::unique_ptr<T_Factory_AbstractProduct> Build ()
                {
                    return std::unique_ptr<T_Factory_AbstractProduct> (new T_ConcreteProduct ());
                }

            private:
                Proxy (const Proxy& proxy) = delete;
                Proxy (Proxy&& proxy) = delete;
                Proxy& operator= (const Proxy& proxy) = delete;
        };


    template <class T_Factory_, class T_ConcreteProduct_>
        Proxy<T_Factory_, T_ConcreteProduct_>::Proxy (const T_Factory_Identifier& identifier)
        {
            // get the factory. First time creates it.
            T_Factory_& factory (T_Factory_::Instance ());
            
            // Insert the builder into the factory
            factory.add (identifier, Proxy <T_Factory_, T_ConcreteProduct_>::Build);
        }
};

#endif
