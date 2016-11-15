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

#ifndef SRC_OBJECTIVE_FUNCTIONALS_OBJECTIVEFUNCTIONAL_H_INCLUDE_GUARD
#define SRC_OBJECTIVE_FUNCTIONALS_OBJECTIVEFUNCTIONAL_H_INCLUDE_GUARD

#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshFunction.h>
#include <dolfin/log/dolfin_log.h>
#include <dolfin/fem/assemble.h>
#include <string>
#include <memory>
#include <dcp/problems/SubdomainType.h>
#include <dcp/expressions/VariableExpression.h>
#include <dcp/objective_functionals/GenericObjectiveFunctional.h>

namespace dcp
{
    /*! \class ObjectiveFunctional ObjectiveFunctional.h
     *  \brief Class for objective functionals.
     *
     *  This class represents a generic functional and contains both its representation and its gradient. It derives
     *  from \c GenericObjectiveFunctional and extends its functionalities to a concrete objective functional.
     *  The class is template-ized over the type of the functional, which must be derived from \c dolfin::Form (and will
     *  typically defined in a ufl file with the keyword \c forms).
     *  
     *  Template arguments are:
     *  \arg T_FunctionalForm_ the functional form type
     */

    template <class T_FunctionalForm_>
        class ObjectiveFunctional : public dcp::GenericObjectiveFunctional
    {
        // ---------------------------------------------------------------------------------------------//  

        public:
            /******************* TYPEDEFS *******************/
            typedef T_FunctionalForm_ T_FunctionalForm;

            /******************* CONSTRUCTORS *******************/
            //! Default constructor is deleted. The class is not default constructable.
            ObjectiveFunctional () = delete;

            //!  Constructor [1]
            /*!
             *  \param mesh the mesh over which the functional is defined
             *  \param gradient the gradient of the functional
             *  
             *  The stored mesh's ownership will be shared between the object and the input argument.
             */
            ObjectiveFunctional (const std::shared_ptr<const dolfin::Mesh> mesh, 
                                 const std::shared_ptr<dcp::GenericExpression>& gradient);

            //!  Constructor [2]
            /*!
             *  \param mesh the mesh over which the functional is defined
             *  \param gradient the gradient of the functional
             *  \param functional the functional itself, which will be used to initialize the protected member
             *  of the class
             *  
             *  The stored mesh's ownership will be shared between the object and the input argument.
             */
            ObjectiveFunctional (const std::shared_ptr<const dolfin::Mesh> mesh, 
                                 const std::shared_ptr<dcp::GenericExpression>& gradient, 
                                 const T_FunctionalForm& functional);

            //! Copy constructor (default shallow copy constructro)
            ObjectiveFunctional (const ObjectiveFunctional<T_FunctionalForm>& rhs) = default;

            
            /******************* DESTRUCTOR *******************/
            //! Destructor
            /*! 
             *  Default destructor, since members of the class are trivially 
             *  destructible.
             */
            virtual ~ObjectiveFunctional () {};


            /******************* GETTERS *******************/
            //! Get reference to the functional [1]. Overrides method in \c GenericObjectiveFunctional
            /*! 
             *  \return a reference to the functional form
             */
            virtual const T_FunctionalForm& functional () const override;

            //! Get reference to the functional gradient [1]. Overrides method in \c GenericObjectiveFunctional
            /*! 
             *  \return a reference to the functional gradient
             */
            virtual const dcp::GenericExpression& gradient () const override;


            //! Get reference to the functional [2]. Overrides method in \c GenericObjectiveFunctional
            /*! 
             *  \return a reference to the functional form
             */
            T_FunctionalForm& functional () override;

            //! Get reference to the functional gradient [2]. Overrides method in \c GenericObjectiveFunctional
            /*! 
             *  \return a reference to the functional gradient
             */
            dcp::GenericExpression& gradient () override;


            /******************* SETTERS *******************/
            //! Set coefficients for the protected member variables. Overrides method in \c GenericObjectiveFunctional
            /*!
             *  Possible values for \c coefficientType are:
             *  \li \c functional to set coefficient of the functional
             *  \li \c gradient to set coefficient of the gradient
             *  
             *  See \c GenericObjectiveFunctional documentation for more details.
             */
            virtual void setCoefficient (const std::string& coefficientType, 
                                         const std::shared_ptr<const dolfin::GenericFunction> coefficientValue,
                                         const std::string& coefficientName) override;

            //! Set integration subdomains for the protected member variable \c functional_. 
            //! Overrides method in \c GenericObjectiveFunctional
            /*!
             *  See \c GenericObjectiveFunctional documentation for more details.
             */
            virtual void setIntegrationSubdomain (std::shared_ptr<const dolfin::MeshFunction<std::size_t>> meshFunction,
                                                   const dcp::SubdomainType& subdomainType) override;


            /******************* METHODS *******************/
            //! Evaluate the stored functional. Overrides method in \c GenericObjectiveFunctional
            /*!
             *  See \c GenericObjectiveFunctional documentation for more details.
             */
            virtual double evaluateFunctional () const override;

            //! Evaluate at given point in given cell. Overrides method in 
            //! \c GenericObjectiveFunctional
            /*!
             *  See \c GenericObjectiveFunctional documentation for more details.
             */
            virtual void evaluateGradient (dolfin::Array<double>& values, 
                                           const dolfin::Array<double>& x, 
                                           const ufc::cell& cell) const override;

            //! Evaluate the stored functional gradient at given point. Overrides method in 
            //! \c GenericObjectiveFunctional
            /*!
             *  See \c GenericObjectiveFunctional documentation for more details.
             */
            virtual void evaluateGradient (dolfin::Array<double>& values, 
                                           const dolfin::Array<double>& x) const override;

            // ---------------------------------------------------------------------------------------------//

        protected:
            //! The functional itself
            T_FunctionalForm functional_;

            //! The gradient of the functional
            std::shared_ptr<dcp::GenericExpression> gradient_;

            // ---------------------------------------------------------------------------------------------//

        private:
    };




    // ============================================================================================== //
    // ==================================== IMPLEMENTATION ========================================== //
    // ============================================================================================== //


    /******************* CONSTRUCTORS *******************/
    template <class T_FunctionalForm>
        ObjectiveFunctional<T_FunctionalForm>::
        ObjectiveFunctional (const std::shared_ptr <const dolfin::Mesh> mesh, 
                             const std::shared_ptr<dcp::GenericExpression>& gradient) :
            GenericObjectiveFunctional (mesh),
            functional_ (mesh),
            gradient_ (gradient) 
    {
        dolfin::log (dolfin::DBG, "ObjectiveFunctional object created");
    }



    template <class T_FunctionalForm>
        ObjectiveFunctional<T_FunctionalForm>::
        ObjectiveFunctional (const std::shared_ptr <const dolfin::Mesh> mesh, 
                             const std::shared_ptr<dcp::GenericExpression>& gradient,
                             const T_FunctionalForm& functional) :
            GenericObjectiveFunctional (mesh),
            functional_ (functional),
            gradient_ (gradient) 
    {
        dolfin::log (dolfin::DBG, "ObjectiveFunctional object created");
    }



    /******************* GETTERS *******************/
    template <class T_FunctionalForm>
        const T_FunctionalForm& ObjectiveFunctional<T_FunctionalForm>::
        functional () const
        {
            return functional_;
        }



    template <class T_FunctionalForm>
        const dcp::GenericExpression& ObjectiveFunctional<T_FunctionalForm>::
        gradient () const
        {
            return *gradient_;
        }



    template <class T_FunctionalForm>
        T_FunctionalForm& ObjectiveFunctional<T_FunctionalForm>::
        functional ()
        {
            return functional_;
        }



    template <class T_FunctionalForm>
        dcp::GenericExpression& ObjectiveFunctional<T_FunctionalForm>::
        gradient ()
        {
            return *gradient_;
        }



    /******************* SETTERS *******************/
    template <class T_FunctionalForm>
        void ObjectiveFunctional<T_FunctionalForm>::
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
                gradient_ -> setCoefficient (coefficientName, coefficientValue);
                dolfin::end ();
            }
            else
            {
                dolfin::warning ("Cannot set coefficient in linear differential problem. Form type \"%s\" unknown",
                                 coefficientType.c_str ());
            }
        }



    template <class T_FunctionalForm>
        void ObjectiveFunctional<T_FunctionalForm>::
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
    template <class T_FunctionalForm>
        double ObjectiveFunctional<T_FunctionalForm>::
        evaluateFunctional () const
        {
            return dolfin::assemble (functional_);
        }



    template <class T_FunctionalForm>
        void ObjectiveFunctional<T_FunctionalForm>::
        evaluateGradient (dolfin::Array<double>& values, 
                          const dolfin::Array<double>& x, 
                          const ufc::cell& cell) const
        {
            gradient_ -> eval (values, x, cell);
        }



    template <class T_FunctionalForm>
        void ObjectiveFunctional<T_FunctionalForm>::
        evaluateGradient (dolfin::Array<double>& values, 
                          const dolfin::Array<double>& x) const
        {
            gradient_ -> eval (values, x);
        }
}

#endif

