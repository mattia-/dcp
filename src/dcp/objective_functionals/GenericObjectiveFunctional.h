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

#ifndef SRC_OBJECTIVE_FUNCTIONALS_GENERICOBJECTIVEFUNCTIONAL_H_INCLUDE_GUARD
#define SRC_OBJECTIVE_FUNCTIONALS_GENERICOBJECTIVEFUNCTIONAL_H_INCLUDE_GUARD

#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshFunction.h>
#include <dolfin/common/Array.h>
#include <dolfin/fem/Form.h>
#include <string>
#include <dcp/expressions/VariableExpression.h>
#include <dcp/problems/SubdomainType.h>

namespace dcp
{
    /*! \class GenericObjectiveFunctional GenericObjectiveFunctional.h
     *  \brief Generic base class for objective functionals.
     *
     *  This class contains the basic interface for a generic objective functional. The class is pure virtual 
     *  and it is intended to be use in order to apply polymorphism for concrete instances of the derived 
     *  \c ObjectiveFunctional class, which is derived for this one
     */

    class GenericObjectiveFunctional
    {
        // ---------------------------------------------------------------------------------------------//  

        public:
            /******************* CONSTRUCTORS *******************/
            //! Default constructor is deleted
            GenericObjectiveFunctional () = delete;

            //! Constructor
            GenericObjectiveFunctional (const std::shared_ptr <const dolfin::Mesh> mesh);
            

            /******************* DESTRUCTOR *******************/

            //! Default destructor
            virtual ~GenericObjectiveFunctional () {};

            
            /******************* GETTERS *******************/
            //! Get const reference to the mesh
            /*! 
             *  \return a const reference to the mesh 
             */
            virtual const dolfin::Mesh& mesh () const;

            //! Get reference to the functional [1]
            /*! 
             *  \return a reference to the functional form
             */
            virtual const dolfin::Form& functional () const = 0;

            //! Get reference to the functional gradient [1]
            /*! 
             *  \return a reference to the functional gradient
             */
            virtual const dolfin::Expression& gradient () const = 0;

            //! Get reference to the functional [2]
            /*! 
             *  \return a reference to the functional form
             */
            virtual dolfin::Form& functional () = 0;

            //! Get reference to the functional gradient [2]
            /*! 
             *  \return a reference to the functional gradient
             */
            virtual dolfin::Expression& gradient () = 0;


            /******************* SETTERS *******************/
            //! Set coefficients for the protected member variables
            /*!
             *  This method is meant to be overridden in derived classes.
             *  Parameters are:
             *  \param coefficientType used to disambiguate between different member variables to choose which coefficient
             *  to set
             *  \param coefficientValue value of the coefficient
             *  \param coefficientName string identifying the coefficient to set
             */
            virtual void setCoefficient (const std::string& coefficientType, 
                                         const std::shared_ptr<const dolfin::GenericFunction> coefficientValue,
                                         const std::string& coefficientName) = 0;

            //! Set integration subdomains for the protected member variable that represent the functional (which must
            //! be declared in the derived class)
            /*! 
             *  This method is meant to be overridden in derived classes.
             *  Input arguments are:
             *  \param meshFunction the mesh function used to set the integration subdomains
             *  \param subdomainType the type of the subdomains, chosen among those provided by the enumeration
             *  class \c dcp::SubdomainType
             */
            virtual void setIntegrationSubdomain (std::shared_ptr<const dolfin::MeshFunction<std::size_t>> meshFunction,
                                                   const dcp::SubdomainType& subdomainType) = 0;

            /******************* METHODS *******************/
            //! Evaluate the stored functional
            /*!
             *  \return a double containing the evaluation of the functional on the mesh
             */
            virtual double evaluateFunctional () const = 0;

            //! Evaluate the stored functional gradient at given point
            /*!
             *  Input arguments are:
             *  \param values the values at the point
             *  \param x the coordinates of the point
             *  \param cell the cell which contains the given point
             *  
             *  \return a double containing the evaluation of the functional on the mesh
             */
            virtual void evaluateGradient (dolfin::Array<double>& values, 
                                           const dolfin::Array<double>& x, 
                                           const ufc::cell& cell) const = 0;

            //! Evaluate at given point in given cell
            /*!
             *  Input arguments are:
             *  \param values the values at the point
             *  \param x the coordinates of the point
             *  
             *  \return a double containing the evaluation of the functional on the mesh
             */
            virtual void evaluateGradient (dolfin::Array<double>& values, 
                                           const dolfin::Array<double>& x) const = 0;

            
            // ---------------------------------------------------------------------------------------------//

        protected:
            //! The mesh over which the functional is defined
            std::shared_ptr <const dolfin::Mesh> mesh_;
            
            
            // ---------------------------------------------------------------------------------------------//

        private:
    };
}
#endif

