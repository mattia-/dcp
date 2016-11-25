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

#ifndef SRC_OPTIMIZERS_DIRICHLETCONTROLUPDATER_H_INCLUDE_GUARD
#define SRC_OPTIMIZERS_DIRICHLETCONTROLUPDATER_H_INCLUDE_GUARD

#include <dcp/problems/EquationSystem.h>
#include <dolfin/function/GenericFunction.h>
#include <dolfin/function/FunctionSpace.h>
#include <string>
#include <dcp/functions/TimeDependentFunction.h>
#include <dcp/subdomains/Subdomain.h>

namespace dcp
{
    /*! \class DirichletControlUpdater DirichletControlUpdater.h
     *  \brief Class to update the value of the control variable in Dirichlet boundary control problems.
     *  
     *  This class is a functor which can be passed to the constructor of any class of the \c GenericImplementer
     *  hierarchy, which will use it to update the primal problem using the new value of the control parameter as the
     *  optimization proceeds.  Note that only the primal system is updated. If the adjoint system needs to be updated
     *  too a new class should be defined (maybe deriving from this one)
     */
    class DirichletControlUpdater
    {
        // ---------------------------------------------------------------------------------------------//

        public:
            /************************* CONSTRUCTORS ********************/
            //! Default constructor is deleted
            DirichletControlUpdater () = delete;
            
            //! Constructor
            /*! 
             *  Input arguments are:
             *  \param problemName string that identifies the problem (in the \c EquationSystem object 
             *  passed as input to <tt> this->operator() ()</tt> ) which contains the control parameter to be updated
             *  \param dirichletBCName the name under which the control Dirichlet boundary condition is stored in the 
             *  protected member \c map in the problem
             *  \param dirichletBoundary the boundary over which the control Dirichlet condition should be enforced
             *  \param component the component of function space on which the dirichlet condition should be enforced;
             *  this allows to enforce boundary conditions on subspaces. A negative value is a placeholder that stands
             *  for a boundary condition on the whole space. A number greater than or equal to 0 will result in a
             *  boundary condition on the subspace with that index number. Default value: -1
             */
            DirichletControlUpdater (const std::string& problemName, 
                                     const std::string& dirichletBCName,
                                     const dcp::Subdomain& dirichletBoundary,
                                     const int& component);


            /************************* OPERATORS ********************/
            //! Call operator [1]
            /*! 
             *  This will actually perform the updating of the primal problem (which is supposed to be the first element
             *  in \c systems )
             *  Input parametes are:
             *  \param systems the systems on which to operate
             *  \param dirichletBCValue the new value for the control Dirichlet boundary condition 
             */
            void operator() (const std::vector<std::shared_ptr<dcp::GenericEquationSystem>> systems,
                             const dolfin::GenericFunction& dirichletBCValue) const;

            //! Call operator [2]
            /*! 
             *  This will actually perform the updating of the primal problem (which is supposed to be the first element
             *  in \c systems )
             *  Input parametes are:
             *  \param systems the systems on which to operate
             *  \param dirichletBCValue the new value for the control Dirichlet boundary condition 
             */
            void operator() (const std::vector<std::shared_ptr<dcp::GenericEquationSystem>> systems,
                             const dcp::TimeDependentFunction& dirichletBCValue) const;

            // ---------------------------------------------------------------------------------------------//

        protected:
            //! The name identifying the problem inside the \c EquationSystem which contains the
            //! control parameter
            std::string problemName_;

            //! The name under which the control Dirichlet boundary condition is stored in the 
            //! protected member \c map in the problem
            std::string dirichletBCName_;

            //! The boundary over which the control Dirichlet condition should be enforced
            const dcp::Subdomain& dirichletBoundary_;

            //! The component of the function space on which the dirichlet condition should be enforced.
            int component_;
    };
}

#endif
