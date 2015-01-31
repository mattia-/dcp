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

#include <dcp/splitting_methods/AbstractSplittingMethod.h>
#include <dcp/differential_problems/TimeDependentEquationSystem.h>


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
    class NavierStokesSplittingMethod : public AbstractSplittingMethod
    {
        // ---------------------------------------------------------------------------------------------//

        public:
            /************************* CONSTRUCTORS ********************/
            //!  Constructor 
            /*!
             *  \param functionSpaces the function spaces of the various differential problems that will be stored in 
             *  the protected member \c differentialSystem_. The first element in the initializer list passed should be
             *  the velocity function space, the second one the pressure function space
             */
            NavierStokesSplittingMethod (std::initializer_list<dolfin::FunctionSpace> functionSpaces);

            
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

            
            /********************** GETTERS ***********************/
            //! Get the system stored in the protected member \c differentialSystem_ [1], const version
            /*! 
             *  \return a pointer to the system
             */
            virtual const dcp::TimeDependentEquationSystem& system () const;

            //! Get the system stored in the protected member \c differentialSystem_ [2], non const version
            /*! 
             *  \return a pointer to the system
             */
            virtual dcp::TimeDependentEquationSystem& system ();

            //! Access problem in \c differentialSystem_ with given name [1] (read only)
            /*!
             *  \param name name of the problem to be accessed. 
             *  
             *  \return a reference to the problem
             */
            virtual const dcp::TimeDependentProblem& problem (const std::string& name) const;
            
            //! Access problem in \c differentialSystem_ with given name [2] (read and write)
            /*!
             *  \param name name of the problem to be accessed. 
             *  
             *  \return a reference to the problem
             */
            virtual dcp::TimeDependentProblem& problem (const std::string& name);
            
            
            // ---------------------------------------------------------------------------------------------//

        protected:
            //! The velocity function space. It is just a pointer to the first element of \c functionSpaces_, which
            //! is inherited from \c dcP::AbstractSplittingMethod
            const std::shared_ptr<dolfin::FunctionSpace> velocityFunctionSpace_;
            
            //! The pressure function space. It is just a pointer to the second element of \c functionSpaces_, which
            //! is inherited from \c dcP::AbstractSplittingMethod
            const std::shared_ptr<dolfin::FunctionSpace> pressureFunctionSpace_;
            

            // ---------------------------------------------------------------------------------------------//

        private:

    };

}
#endif
