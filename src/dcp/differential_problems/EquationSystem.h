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
            virtual void solve () override;
            
            //! Solve the problem corresponding to the name given [1]
            /*!
             *  \param problemName a string identifying the problem to be solved. If no problem with that name
             *  is found, a warning is issued
             */
            virtual void solve (const std::string& problemName) override;
            
        // ---------------------------------------------------------------------------------------------//  

        protected:
    };
}

#endif
