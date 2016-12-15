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

#ifndef SRC_OPTIMIZERS_GENERICDESCENTMETHOD_H_INCLUDE_GUARD
#define SRC_OPTIMIZERS_GENERICDESCENTMETHOD_H_INCLUDE_GUARD

#include <dolfin/parameter/Parameters.h>
#include <dolfin/function/GenericFunction.h>
#include <dolfin/function/Function.h>
#include <dolfin/function/Expression.h>
#include <dcp/objective_functionals/GenericObjectiveFunctional.h>
#include <dcp/problems/GenericEquationSystem.h>
#include <dcp/optimizers/GradientSearchDirection.h>
#include <dcp/utils/DotProduct.h>
#include <functional>

namespace dcp
{
    /*! \class GenericDescentMethod GenericDescentMethod.h
     *  \brief Generic base class for descent methods.
     *
     *  This class defines the base interface for all descent methods.
     *  It provides a \c apply() method to perform the optimization of the
     *  problem passed as input argument and an empty parameters set that
     *  can be populated by derived classes to store concrete methods' settings
     */

    class GenericDescentMethod
    {
        // ---------------------------------------------------------------------------------------------//

        public:
            /************************* CONSTRUCTORS ********************/
            //! Default constructor
            GenericDescentMethod ();


            /************************* DESTRUCTOR ********************/
            //! Destructor
            /*! Default destructor, since members of the class are trivially
             *  destructible.
             *  The protected member will all be initialized with default parameters (see specific members documentation
             *  for details).
             */
            virtual ~GenericDescentMethod () {};

            /********************** VARIABLES ***********************/
            //! the problem parameters
            dolfin::Parameters parameters;

            // ---------------------------------------------------------------------------------------------//

        protected:

            // ---------------------------------------------------------------------------------------------//

        private:

    };
}

#endif
