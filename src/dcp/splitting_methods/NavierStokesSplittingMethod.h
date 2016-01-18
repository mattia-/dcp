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

#ifndef SRC_SPLITTING_METHODS_NAVIERSTOKESSPLITTINGMETHOD_H_INCLUDE_GUARD
#define SRC_SPLITTING_METHODS_NAVIERSTOKESSPLITTINGMETHOD_H_INCLUDE_GUARD

#include <dcp/splitting_methods/GenericSplittingMethod.h>
#include <dcp/problems/TimeDependentEquationSystem.h>


namespace dcp
{
    /*! \class NavierStokesSplittingMethod NavierStokesSplittingMethod.h
     *  \brief Base class for splitting methods for a Navier-Stokes problem.
     *         
     *  This class contains the basic interface for a splitting method to 
     *  solve a Navier-Stokes problem
     *  It is an abstract class, it only provides the basic interface and needs
     *  to be derived to implement concrete splitting methods, like Chorin-Temam
     *  and Incremental Chorin-Temam. 
     */
    class NavierStokesSplittingMethod : public GenericSplittingMethod
    {
        // ---------------------------------------------------------------------------------------------//

        public:
            /************************* CONSTRUCTORS ********************/
            //!  Constructor 
            /*!
             *  \param functionSpaces the function spaces of the various differential problems that will be stored in 
             *  the protected member \c differentialSystem_. The first element in the vector passed should be a pointer
             *  to the velocity function space, the second one a pointer to the pressure function space
             */
            NavierStokesSplittingMethod (const std::vector<std::shared_ptr <dolfin::FunctionSpace>> functionSpaces);
            
            /************************* DESTRUCTOR ********************/
            //! Destructor
            /*! 
             *  Default destructor, since members of the class are trivially destructible. 
             *  It is declared as pure virtual to make the class abstract. Indeed, one cannot use this class because
             *  the protected member \c differentialSystem_ is not set (or better, it is set to \c nullptr and thus
             *  using the class methods would result in a segmentation fault). On the other hand, all methods can
             *  be defined here and inherited in derived classes. That's why we defined a pure virtual destructor.
             *  Note that it still needs to be defined, as it will be called anyway when the base class is destroyed.
             */
            virtual ~NavierStokesSplittingMethod () = 0;

            // ---------------------------------------------------------------------------------------------//

        protected:
            //! The velocity function space. It is just a reference to the object pointed by the first element of 
            //! \c functionSpaces_, which is inherited from \c dcp::GenericSplittingMethod
            const dolfin::FunctionSpace& velocityFunctionSpace_;
            
            //! The pressure function space. It is just a reference to the object pointed by the second element of 
            //! \c functionSpaces_, which is inherited from \c dcp::GenericSplittingMethod
            const dolfin::FunctionSpace& pressureFunctionSpace_;
            

            // ---------------------------------------------------------------------------------------------//

        private:

    };

}
#endif
