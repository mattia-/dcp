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

#ifndef SRC_SUBDOMAINS_SUBDOMAIN_H_INCLUDE_GUARD
#define SRC_SUBDOMAINS_SUBDOMAIN_H_INCLUDE_GUARD

#include <dolfin/mesh/SubDomain.h>
#include <functional>

namespace dcp
{
    /*! class Subdomain Subdomain.h
     *  \brief Class for subdomains definition.
     *  
     *  This class offers an alernative to the dolfin way to define subdomains. 
     *  Instead of having the users derive a class from \c dolfin::SubDomain and override the \c inside() method,
     *  this class stores a function wrapper as protected member, initialized upon building through the 
     *  object passed to the constructor. This function wrapper is called when the method \c inside() is called to
     *  decide wether a point is inside the subdomain or not.
     */
    class Subdomain : public dolfin::SubDomain
    {
        // ---------------------------------------------------------------------------------------------//  
        public:
            typedef std::function <bool (const dolfin::Array<double>&, bool)> Evaluator;
            
            /******************* CONSTRUCTORS *******************/
            //! Default constructor
            /*
             *  Input arguments
             *  \param evaluator the evaluator to be used when calling the \c inside() method
             */
            Subdomain (const Evaluator& evaluator);
            
            
            /******************* DESTRUCTOR *******************/
            //! Default destructor                
            /*! 
             *  Default destructor, since members of the class are trivially 
             *  destructible.
             */
            virtual ~Subdomain () {};


            /******************* METHODS *******************/
            //! Return \c true for points inside the domain
            /*!
             *  Input arguments are:
             *  \param x the coordinates of the point
             *  \param on_boundary \c true for points on the boundary
             */
            virtual bool inside (const dolfin::Array<double>& x, bool on_boundary) const;

        // ---------------------------------------------------------------------------------------------//  
        protected:
            //! The evaluator to use when the \c inside() method is called
            Evaluator evaluator_;
    };
}

#endif
