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

#ifndef SRC_OPTIMIZERS_EMPTYUPDATER_H_INCLUDE_GUARD
#define SRC_OPTIMIZERS_EMPTYUPDATER_H_INCLUDE_GUARD

#include <dcp/problems/EquationSystem.h>
#include <dolfin/function/GenericFunction.h>
#include <dolfin/function/FunctionSpace.h>
#include <string>
#include <dcp/functions/TimeDependentFunction.h>
#include <dcp/subdomains/Subdomain.h>

namespace dcp
{
    /*! \class EmptyUpdater EmptyUpdater.h
     *  \brief Class to update control problems which actually does nothing.
     *  
     *  This class is a functor which can be passed to the constructor of any class of the \c GenericImplementer
     *  hierarchy, which will use it to update the primal problem using the new value of the control parameter as the
     *  optimization proceeds.  Note that only the primal system is updated. If the adjoint system needs to be updated
     *  too a new class should be defined (maybe deriving from this one). 
     *  This class actually does nothing, but it is provided for compatibility in those cases where no action is needed
     *  for the update of the control problems.
     */
    class EmptyUpdater
    {
        // ---------------------------------------------------------------------------------------------//

        public:
            /************************* CONSTRUCTORS ********************/
            //! Default constructor 
            EmptyUpdater () = default;

            /************************* OPERATORS ********************/
            //! Call operator [1]
            /*! 
             *  This will do nothing
             *
             *  Input parametes are:
             *  \param systems the systems on which to operate
             *  \param dirichletBCValue the new value for the control Dirichlet boundary condition 
             */
            void operator() (const std::vector<std::shared_ptr<dcp::GenericEquationSystem>> systems,
                             const dolfin::GenericFunction& dirichletBCValue) const;

            //! Call operator [2]
            /*! 
             *  This will do nothing
             *
             *  Input parametes are:
             *  \param systems the systems on which to operate
             *  \param dirichletBCValue the new value for the control Dirichlet boundary condition 
             */
            void operator() (const std::vector<std::shared_ptr<dcp::GenericEquationSystem>> systems,
                             const dcp::TimeDependentFunction& dirichletBCValue) const;

            // ---------------------------------------------------------------------------------------------//

        protected:
    };
}

#endif

