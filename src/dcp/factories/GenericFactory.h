/* 
 *  Copyright (C) 2014, Mattia Tamellini, mattia.tamellini@gmail.com
 * 
 *  This file is part of the DCP library
 *   
 *   The DCP library is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   The DCP library is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with the DCP library.  If not, see <http://www.gnu.org/licenses/>. 
 */ 

#ifndef SRC_FACTORY_GENERICFACTORY_H_INCLUDE_GUARD
#define SRC_FACTORY_GENERICFACTORY_H_INCLUDE_GUARD

#include <memory>
#include <functional>
#include <map>
#include <vector>
#include <iostream>
#include <string>

namespace dcp
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
            class T_GenericProduct_, 
            class T_Identifier_, 
            class T_Builder_ = std::function <std::unique_ptr<T_GenericProduct_> ()> 
        >
            class GenericFactory
            {
                public:
                    typedef T_GenericProduct_ T_GenericProduct;
                    typedef T_Identifier_      T_Identifier;
                    typedef T_Builder_         T_Builder;

                    //! Method to access the only instance of the factory
                    static GenericFactory& Instance();

                    //! Get the object with given identifier 
                    /*!
                      The pointer is null if no match was found for the identifier.
                      */
                    std::unique_ptr<T_GenericProduct> create (const T_Identifier& identifier) const;

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
            class T_GenericProduct,
            class T_Identifier,
            class T_Builder
        >
        GenericFactory<T_GenericProduct, T_Identifier, T_Builder>& 
        GenericFactory<T_GenericProduct, T_Identifier, T_Builder>::Instance() 
        {
            static GenericFactory factory;
            return factory;
        }
    


    template
        <
            class T_GenericProduct,
            class T_Identifier,
            class T_Builder
        >
        std::unique_ptr <T_GenericProduct> 
        GenericFactory<T_GenericProduct, T_Identifier, T_Builder>::create(const T_Identifier& identifier) const
        {
            auto f = storedData_.find (identifier); 
            if (f == storedData_.end ())
            {
                dolfin::warning ("identifier %s not found in factory", identifier.c_str ());
                return nullptr;
            }
            else
            {
                return std::unique_ptr<T_GenericProduct> (f->second());
            }
        }
    


    template
        <
            class T_GenericProduct,
            class T_Identifier,
            class T_Builder
        >
        void 
        GenericFactory<T_GenericProduct, T_Identifier, T_Builder>::
        add (const T_Identifier& identifier, const T_Builder& builder)
        {
            auto f = storedData_.insert(std::make_pair(identifier, builder));
            if (f.second == false)
            {
                dolfin::warning ("double registration in factory while trying to register product %s", identifier.c_str ());
            }
        }


    
    template
        <
            typename T_GenericProduct,
            typename T_Identifier,
            typename T_Builder
        >
        std::vector<T_Identifier> 
        GenericFactory<T_GenericProduct, T_Identifier, T_Builder>::registered () const
        {
            std::vector<T_Identifier> tmp;
            tmp.reserve (storedData_.size());
            
            for(auto i = storedData_.begin(); i != storedData_.end(); ++i)
                tmp.push_back (i->first);
            
            return tmp;
        }
};

#endif
