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

#ifndef SRC_OPTIMIZERS_GRADIENTSEARCHDIRECTION_H_INCLUDE_GUARD
#define SRC_OPTIMIZERS_GRADIENTSEARCHDIRECTION_H_INCLUDE_GUARD

#include <dolfin/function/GenericFunction.h>
#include <dolfin/function/Function.h>

namespace dcp
{
    /*! \class GradientSearchDirection GradientSearchDirection.h
     *  \brief Helper class that defines the search direction for descent methods as the opposite of the gradient.
     * 
     *  This class defines an object whose call operator computes the search direction for descent methods by
     *  taking the opposite of the gradient.
     *  It is used as default value for \c searchDirectionComputer_ in \c dcp::BacktrackingOptimizer
     */
    
    class GradientSearchDirection
    {
        // ---------------------------------------------------------------------------------------------//
        
        public:
            /************************* CONSTRUCTORS ********************/
            //! Default constructor
            GradientSearchDirection () = default;
            

            /************************* DESTRUCTOR ********************/
            //! Destructor
            /*! Default destructor, since members of the class are trivially 
             * destructible.
             */
            ~GradientSearchDirection () {};
            
            
            /********************** METHODS ***********************/
            //! Call operator that computes the search direction
            /*!
             *  \param searchDirection will contain the search direction at the end of the function
             *  \param gradient the current value of the gradient of the functional
             */
            void operator() (dolfin::Function& searchDirection, const dolfin::Function gradient);
            
            // ---------------------------------------------------------------------------------------------//

        protected:

            // ---------------------------------------------------------------------------------------------//

        private:
            
    };
}

#endif


