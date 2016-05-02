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

#ifndef SRC_UTILS_DOTPRODUCT_H_INCLUDE_GUARD
#define SRC_UTILS_DOTPRODUCT_H_INCLUDE_GUARD

#include <dolfin/log/dolfin_log.h>
#include <dolfin/function/Function.h>
#include <dolfin/fem/assemble.h>
#include <dolfin/fem/Form.h>
#include <functional>
#include <string>

namespace dcp
{
    /*! \class DotProduct DotProduct.h
     *  \brief Class to compute dot products between finite element functions
     *  
     *  This class allows to compute the dot product of two finite element functions and the norm of a single finite 
     *  element function (by using the well-known relation between dot products and norms). It tries to automatically 
     *  determine the right dot product to use, but it can be overridden in particular cases where such automatic
     *  determination fails. All the implemented forms use piecewise-linear finite element spaces.
     */
    class DotProduct
    {
        // ---------------------------------------------------------------------------------------------//
        
        public:
            /************************* CONSTRUCTORS ********************/
            //! Default constructor
            DotProduct ();
            
            /************************* DESTRUCTOR ********************/
            //! Destructor
            /*! Default destructor, since members of the class are trivially 
             * destructible.
             */
            virtual ~DotProduct () {};
            
            
            /********************** METHODS ***********************/
            //! Set the form used to compute dot products and norms. 
            /*!
             *  \param dotProductComputer the form to be stored in the protected member \c dotProductComputer_ and
             *  to be used when \c compute() is called
             */
            virtual void setDotProductComputer (const dolfin::Form& dotProductComputer);
            
            //! Reset the value of the protected membet \c dotProductComputer_ so that the default form will be used
            virtual void resetDotProductComputer ();
                
            //! Function to compute the dot product between two \c dolfin::GenericFunction objects. 
            /*! 
             *  The \c dolfin::Form to be used is the one stored in dotProductComputer_ , if such variable is not empty. 
             *  Otherwise, it will be determined based on the mesh (third input argument). If no mesh is provided, the
             *  mesh of the first input argument is used
             *  \param first the first function of the dot product
             *  \param second the second function of the dot product
             *  \param mesh the mesh on which to perform the dot product
             *  
             *  \return the value of the dot product
             */
            virtual double compute (const dolfin::GenericFunction& first, 
                                    const dolfin::GenericFunction& second,
                                    const dolfin::Mesh& mesh);
            
            //! Function to compute the dot product between two \c dolfin::Function objects. 
            /*! 
             *  The \c dolfin::Form to be used is the one stored in dotProductComputer_ , if such variable is not empty. 
             *  Otherwise, it will be determined based on the mesh (third input argument). If no mesh is provided, the
             *  mesh of the first input argument is used
             *  \param first the first function of the dot product
             *  \param second the second function of the dot product
             *  
             *  The first function mesh will be used
             *  
             *  \return the value of the dot product
             */
            virtual double compute (const dolfin::Function& first, 
                                    const dolfin::Function& second);
            
            //! Function to compute the norm of a \c dolfin::GenericFunction object.
            /*! 
             *  This function simply wraps a call to \c compute().
             *  \param function the function or expression of which to compute the norm
             *  \param mesh the mesh to be used
             */
            virtual double norm (const dolfin::GenericFunction& function, 
                                 const dolfin::Mesh& mesh);
            
            //! Function to compute the norm of a \c dolfin::Function object.
            /*! 
             *  This function simply wraps a call to \c compute().
             *  \param function the function or expression of which to compute the norm
             */
            virtual double norm (const dolfin::Function& function);
            
            // ---------------------------------------------------------------------------------------------//

        protected:
            /********************** METHODS ***********************/
            //! Function to get the right dotProductComputer
            /*!
             *  \param first the first function of the dot product
             *  \param second the second function of the dot product
             *  \param mesh the mesh on which to perform the dot product 
             *  
             *  \return the dotProductComputer itself
             */
            std::shared_ptr<dolfin::Form> getDotProductComputer_ (const dolfin::GenericFunction& first,
                                                                  const dolfin::GenericFunction& second,
                                                                  const dolfin::Mesh& mesh);
            
            /********************** MEMBERS ***********************/
            //! The form that will be used to compute the dot product
            /*! 
             *  The default value is \c nullptr, in which case the function will use one  of the forms in 
             *  \c dotproductforms.h, trying to determine the right one by checking the geometrical dimensions
             *  of the input objects. However, sometimes it may be useful to have a user-defined object to perform
             *  the task. To do so, use the function \c setDotProductComputer
             */
            std::shared_ptr<dolfin::Form> dotProductComputer_;
            
            // ---------------------------------------------------------------------------------------------//

        private:
            
    };
}

#endif


