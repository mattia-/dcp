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

#ifndef SRC_OBJECTIVE_FUNCTIONAL_OBJECTIVEFUNCTIONAL_H_INCLUDE_GUARD
#define SRC_OBJECTIVE_FUNCTIONAL_OBJECTIVEFUNCTIONAL_H_INCLUDE_GUARD

#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshFunction.h>
#include <dolfin/log/dolfin_log.h>
#include <dolfin/fem/assemble.h>
#include <string>
#include <memory>
#include <dcp/differential_problems/SubdomainType.h>
#include <dcp/expressions/VariableExpression.h>
#include <dcp/objective_functional/AbstractObjectiveFunctional.h>

namespace dcp
{
    /*! \class ObjectiveFunctional ObjectiveFunctional.h
     *  \brief Class for objective functionals.
     *
     *  This class represents a generic functional and contains both its representation and its gradient. It derives
     *  from \c AbstractObjectiveFunctional and extends its functionalities to a concrete objective functional.
     *  The class is template-ized over the type of the functional, which must be derived from \c dolfin::Form (and will
     *  typically defined in a ufl file with the keyword \c forms).
     *  
     *  Template arguments are:
     *  \arg T_FunctionalForm_ the functional form type
     *  \arg T_Gradient_ the gradient type, which must be a class derived from \c dcp::GenericExpression
     */

    template <class T_FunctionalForm_, class T_Gradient_>
        class ObjectiveFunctional : public dcp::AbstractObjectiveFunctional
    {
        // ---------------------------------------------------------------------------------------------//  

        public:
            /******************* TYPEDEFS *******************/
            typedef T_FunctionalForm_ T_FunctionalForm;
            typedef T_Gradient_       T_Gradient;

            /*  
             *  TODO I need another set of constructors that take also a pointer to dcp::Time in input so that we
             *  can build time dependent expressions, which have a different constuctor signature
             *  
             *  Even though for the time being the constructor using the T_Gradient& object can be enough
             *  
             *  Also because then I would need to add a set of constructors for every constructor argument in
             *  dcp::GenericExpression hierarchy...and they are a lot! So maybe I'll be fine with just what I have.
             */

            /******************* CONSTRUCTORS *******************/
            //! Default constructor is deleted. The class is not default constructable.
            ObjectiveFunctional () = delete;

            //!  Constructor from shared pointer [1]
            /*!
             *  \param mesh the mesh over which the functional is defined
             *  \param evaluator the evaluator to be used to build \c gradient_
             *  
             *  The stored mesh's ownership will be shared between the object and the input argument.
             *  The gradient will be created passing \c evaluator to the constructor.
             *  The functional form will be created too, calling the constructor which takes the mesh
             *  as input.
             */
            ObjectiveFunctional (const std::shared_ptr <dolfin::Mesh> mesh
                                 const T_Gradient::Evaluator evaluator);


            //! Constructor from reference [1]
            /*!
             *  \param mesh the mesh over which the functional is defined
             *  \param evaluator the evaluator to be used to build \c gradient_
             *  
             *  The stored mesh's ownership will be unique to the object, since the protected member \c mesh_ is
             *  initialized using the \c new operator and mesh's copy constructor.
             *  The gradient will be created passing \c evaluator to the constructor.
             *  The functional form will be created too, calling the constructor which takes the mesh
             *  as input.
             */
            ObjectiveFunctional (const dolfin::Mesh& mesh,
                                 const T_Gradient::Evaluator evaluator);

            //!  Constructor from shared pointer [2]
            /*!
             *  \param mesh the mesh over which the functional is defined
             *  \param gradient the gradient of the functional
             *  \param functional the functional itself, which will be used to initialize the protected member
             *  of the class
             *  
             *  The stored mesh's ownership will be shared between the object and the input argument.
             */
            ObjectiveFunctional (const std::shared_ptr <dolfin::Mesh> mesh, 
                                 const T_Gradient& gradient, 
                                 const T_FunctionalForm& functional);


            //! Constructor from reference [2]
            /*!
             *  \param mesh the mesh over which the functional is defined
             *  \param gradient the gradient of the functional
             *  \param functional the functional itself, which will be used to initialize the protected member
             *  of the class
             *  
             *  The stored mesh's and ownership will be unique to the object, since the protected member \c mesh_ is
             *  initialized using the \c new operator and mesh's copy constructor.
             */
            ObjectiveFunctional (const dolfin::Mesh& mesh, 
                                 const T_Gradient& gradient, 
                                 const T_FunctionalForm& functional);

            //! Copy constructor
            /*!
             *  Construct from object of type ObjectiveFunctional<T_FunctionalForm, T_Gradient>.
             *  Input arguments are:
             *  \param rhs the object to copy
             *  \param copyMode the type of copy to be performed. 
             *  It can be either \c deep_copy or \c shallow_copy. If the former is selected, the new
             *  object will be created calling the constructor that takes a refenrence to \c dolfin::Mesh, so that the
             *  mesh itself is created. If the latter is selected, the pointer will be copied, adding the new 
             *  ObjectiveFunctional object to the ownership pool of the mesh object. Default value: \c shallow_copy
             */
            ObjectiveFunctional (const ObjectiveFunctional<T_FunctionalForm, T_Gradient>& rhs,
                                 const std::string& copyMode = "shallow_copy");

            
            /******************* DESTRUCTOR *******************/
            //! Destructor
            /*! 
             *  Default destructor, since members of the class are trivially 
             *  destructible.
             */
            virtual ~ObjectiveFunctional () {};


            /******************* GETTERS *******************/
            //! Get const reference to the functional. Overrides method in \c AbstractObjectiveFunctional
            /*! 
             *  \return a const reference to the functional form
             */
            virtual const T_FunctionalForm& functional () const override;

            //! Get const reference to the functional gradient. Overrides method in \c AbstractObjectiveFunctional
            /*! 
             *  \return a const reference to the functional gradient
             */
            virtual const T_Gradient& gradient () const override;


            /******************* SETTERS *******************/
            //! Set coefficients for the protected member variables. Overrides method in \c AbstractObjectiveFunctional
            /*!
             *  Possible values for \c coefficientType are:
             *  \li \c functional to set coefficient of the functional
             *  \li \c gradient to set coefficient of the gradient
             *  
             *  See \c AbstractObjectiveFunctional documentation for more details.
             */
            virtual void setCoefficient (const std::string& coefficientType, 
                                         const std::shared_ptr<const dolfin::GenericFunction> coefficientValue,
                                         const std::string& coefficientName) override;

            //! Set integration subdomains for the protected member variable \c functional_. 
            //! Overrides method in \c AbstractObjectiveFunctional
            /*!
             *  See \c AbstractObjectiveFunctional documentation for more details.
             */
            virtual void setIntegrationSubdomain (std::shared_ptr<const dolfin::MeshFunction<std::size_t>> meshFunction,
                                                   const dcp::SubdomainType& subdomainType) override;


            /******************* METHODS *******************/
            //! Evaluate the stored functional. Overrides method in \c AbstractObjectiveFunctional
            /*!
             *  See \c AbstractObjectiveFunctional documentation for more details.
             */
            virtual double evaluateFunctional () const override;

            //! Evaluate at given point in given cell. Overrides method in 
            //! \c AbstractObjectiveFunctional
            /*!
             *  See \c AbstractObjectiveFunctional documentation for more details.
             */
            virtual void evaluateGradient (dolfin::Array<double>& values, 
                                           const dolfin::Array<double>& x, 
                                           const ufc::cell& cell) const override;

            //! Evaluate the stored functional gradient at given point. Overrides method in 
            //! \c AbstractObjectiveFunctional
            /*!
             *  See \c AbstractObjectiveFunctional documentation for more details.
             */
            virtual void evaluateGradient (dolfin::Array<double>& values, 
                                           const dolfin::Array<double>& x) const override;

            // ---------------------------------------------------------------------------------------------//

        protected:
            //! The functional itself
            T_FunctionalForm functional_;

            //! The gradient of the functional
            T_Gradient gradient_;

            // ---------------------------------------------------------------------------------------------//

        private:
    };




    // ============================================================================================== //
    // ==================================== IMPLEMENTATION ========================================== //
    // ============================================================================================== //


    /******************* CONSTRUCTORS *******************/
    template <class T_FunctionalForm, class T_Gradient>
        ObjectiveFunctional<T_FunctionalForm, T_Gradient>::
        ObjectiveFunctional (const std::shared_ptr <dolfin::Mesh> mesh
                             const T_Gradient::Evaluator evaluator) :
            AbstractObjectiveFunctional (mesh),
            functional_ (*mesh),
            gradient_ (evaluator)
    {
        dolfin::begin (dolfin::DBG, "Creating ObjectiveFunctional...");
        dolfin::log (dolfin::DBG, "ObjectiveFunctional object created");
        dolfin::end ();
    }



    template <class T_FunctionalForm, class T_Gradient>
        ObjectiveFunctional<T_FunctionalForm, T_Gradient>::
        ObjectiveFunctional (const dolfin::Mesh& mesh,
                             const T_Gradient::Evaluator evaluator) :
            AbstractObjectiveFunctional (mesh),
            functional_ (mesh),
            gradient_ (evaluator)
    {
        dolfin::begin (dolfin::DBG, "Creating ObjectiveFunctional...");
        dolfin::log (dolfin::DBG, "ObjectiveFunctional object created");
        dolfin::end ();
    }



    template <class T_FunctionalForm, class T_Gradient>
        ObjectiveFunctional<T_FunctionalForm, T_Gradient>::
        ObjectiveFunctional (const std::shared_ptr <dolfin::Mesh> mesh, 
                             const T_Gradient& gradient, 
                             const T_FunctionalForm& functional) : 
            AbstractObjectiveFunctional (mesh),
            functional_ (functional),
            gradient_ (gradient) // TODO remember it is a shallow copy
    {
        dolfin::begin (dolfin::DBG, "Creating ObjectiveFunctional...");
        dolfin::log (dolfin::DBG, "ObjectiveFunctional object created");
        dolfin::end ();
    }



    template <class T_FunctionalForm, class T_Gradient>
        ObjectiveFunctional<T_FunctionalForm, T_Gradient>::
        ObjectiveFunctional (const dolfin::Mesh& mesh, 
                             const T_Gradient& gradient, 
                             const T_FunctionalForm& functional) : 
            AbstractObjectiveFunctional (mesh),
            functional_ (functional),
            gradient_ (gradient) // TODO remember it is a shallow copy
    {
        dolfin::begin (dolfin::DBG, "Creating ObjectiveFunctional...");
        dolfin::log (dolfin::DBG, "ObjectiveFunctional object created");
        dolfin::end ();
    }



    template <class T_FunctionalForm, class T_Gradient>
        ObjectiveFunctional<T_FunctionalForm, T_Gradient>::
        ObjectiveFunctional (const ObjectiveFunctional<T_FunctionalForm, T_Gradient>& rhs, 
                             const std::string& copyMode) :
            AbstractObjectiveFunctional (copyMode == "shallow_copy" ? 
                                         rhs.mesh_ : 
                                         std::shared_ptr<const dolfin::Mesh> (new dolfin::Mesh (*(rhs.mesh_)))),
            functional_ (rhs.functional_),
            gradient_ (rhs.gradient_) // TODO remember it is a shallow copy
    {
        dolfin::begin (dolfin::DBG, "Creating ObjectiveFunctional ...");
        dolfin::log (dolfin::DBG, "Copy of ObjectiveFunctional created");
        dolfin::end ();
    }



    /******************* GETTERS *******************/
    template <class T_FunctionalForm, class T_Gradient>
        const T_FunctionalForm& ObjectiveFunctional<T_FunctionalForm, T_Gradient>::
        functional () const
        {
            return functional_;
        }



    template <class T_FunctionalForm, class T_Gradient>
        const T_Gradient& ObjectiveFunctional<T_FunctionalForm, T_Gradient>::
        gradient () const
        {
            return gradient_;
        }



    /******************* SETTERS *******************/
    template <class T_FunctionalForm, class T_Gradient>
        void ObjectiveFunctional<T_FunctionalForm, T_Gradient>::
        setCoefficient (const std::string& coefficientType, 
                        const std::shared_ptr<const dolfin::GenericFunction> coefficientValue,
                        const std::string& coefficientName)
        {
            if (coefficientType == "functional")
            {
                dolfin::begin (dolfin::DBG, "Setting functional coefficient \"%s\"...", coefficientName.c_str ());
                functional_.set_coefficient (coefficientName, coefficientValue);
                dolfin::end ();
            }
            else if (coefficientType == "gradient")
            {
                dolfin::begin (dolfin::DBG, "Setting gradient coefficient \"%s\"...", coefficientName.c_str ());
                gradient_.setCoefficient (coefficientName, coefficientValue);
                dolfin::end ();
            }
            else
            {
                dolfin::warning ("Cannot set coefficient in linear differential problem. Form type \"%s\" unknown",
                                 coefficientType.c_str ());
            }
        }



    template <class T_FunctionalForm, class T_Gradient>
        void ObjectiveFunctional<T_FunctionalForm, T_Gradient>::
        setIntegrationSubdomain (std::shared_ptr<const dolfin::MeshFunction<std::size_t>> meshFunction,
                                  const dcp::SubdomainType& subdomainType)
        {
            if (subdomainType == dcp::SubdomainType::INTERNAL_CELLS)
            {
                dolfin::log (dolfin::DBG, "Setting functional integration subdomain on INTERNAL_CELLS...");
                functional_.set_cell_domains (meshFunction);
            }
            else if (subdomainType == dcp::SubdomainType::INTERNAL_FACETS)
            {
                dolfin::log (dolfin::DBG, "Setting functional integration subdomain on INTERNAL_FACETS...");
                functional_.set_interior_facet_domains (meshFunction);
            }
            else if (subdomainType == dcp::SubdomainType::BOUNDARY_FACETS)
            {
                dolfin::log (dolfin::DBG, "Setting functional integration subdomain on BOUNDARY_FACETS...");
                functional_.set_exterior_facet_domains (meshFunction);
            }
            else
            {
                dolfin::warning ("unknown subdomain type requested while trying to apply mesh function to functional");
            }
        }



    /******************* METHODS *******************/
    template <class T_FunctionalForm, class T_Gradient>
        double ObjectiveFunctional<T_FunctionalForm, T_Gradient>::
        evaluateFunctional () const
        {
            return dolfin::assemble (functional_);
        }



    template <class T_FunctionalForm, class T_Gradient>
        void ObjectiveFunctional<T_FunctionalForm, T_Gradient>::
        evaluateGradient (dolfin::Array<double>& values, 
                          const dolfin::Array<double>& x, 
                          const ufc::cell& cell) const
        {
            gradient_.eval (values, x, cell);
        }



    template <class T_FunctionalForm, class T_Gradient>
        void ObjectiveFunctional<T_FunctionalForm, T_Gradient>::
        evaluateGradient (dolfin::Array<double>& values, 
                          const dolfin::Array<double>& x) const
        {
            gradient_.eval (values, x);
        }
}
#endif

