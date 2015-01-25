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

#ifndef SRC_DIFFERENTIAL_PROBLEMS_EQUATIONSYSTEM_H_INCLUDE_GUARD
#define SRC_DIFFERENTIAL_PROBLEMS_EQUATIONSYSTEM_H_INCLUDE_GUARD

#include <dcp/differential_problems/AbstractEquationSystem.h>
#include <map>
#include <tuple>
#include <memory>
#include <iostream>
#include <utility>
#include <vector>
#include <functional>

namespace dcp
{
    /*! \class EquationSystem EquationSystem.h
     *  \brief Class for multi-variable and multi-equation coupled system
     *  
     *  This class derives from AbstractEquationSystem and expands its functionalities by defining the
     *  solve method. This class only works for steady equation systems. For time dependent systems, 
     *  see \c dcp::TimeDependentEquationSystem
     */
    class EquationSystem : public dcp::AbstractEquationSystem
    {
        // ---------------------------------------------------------------------------------------------//  

        public:
            /******************* CONSTRUCTORS ******************/
            //! Default constructor 
            EquationSystem ();
            
            
            /******************* DESTRUCTOR *******************/
            //! Default destructor
            virtual ~EquationSystem () = default;
            

            /******************** METHODS *********************/
            //! Solve all the problems in the order specified by the private member \c solveOrder_
            /*!
             *  \param forceRelinking a boolean flag which, if set to \c true, overrides the current value of protected 
             *  member variable needsLinksScanning_. Default value is \c false
             */
            virtual void solve (const bool& forceRelinking = false) override;
            
            //! Solve the problem corresponding to the name given [1]
            /*!
             *  \param problemName a string identifying the problem to be solved. If no problem with that name
             *  is found, a warning is issued
             *  \param forceRelinking a boolean flag which, if set to \c true, overrides the current value of protected 
             *  member variable needsLinksScanning_. Default value is \c false
             */
            virtual void solve (const std::string& problemName, const bool& forceRelinking = false) override;
            
            //! Solve the problem corresponding to the name given [2]
            /*!
             *  This method is provided only to allow calls like
             *  \code
             *  solve ("foo_problem");
             *  \endcode
             *  In this case, the compiler would in fact otherwise call \c solve \c (const \c bool&) which is the 
             *  best-matching implicit conversion for a parameter of type \c const \c char*. Using this method,
             *  the version of \c solve that takes a \c std::string is called as expected.
             */
            virtual void solve (const char* problemName, const bool& forceRelinking = false) override;
            
        // ---------------------------------------------------------------------------------------------//  

        protected:
    };
}

#endif
