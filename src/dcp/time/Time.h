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

#ifndef SRC_TIMES_TIME_H_INCLUDE_GUARD
#define SRC_TIMES_TIME_H_INCLUDE_GUARD

namespace dcp
{
    /*! class Time Time.h
     *  \brief Class for time definition and management.
     *  
     *  This class defines an object that represents the time during the simulations.
     *  It is not a singleton, so that one may have more than one time object in their program
     *  (even though I'm not sure what good would that be), but the classes in the dcp library
     *  that depend in any way on time all work on objects of type \c dcp::Time& so that when 
     *  using the library a single time object is needed to describe all objects that refer to
     *  the same time value.
     */
    class Time
    {
        // ---------------------------------------------------------------------------------------------//  
        public:
            /******************* CONSTRUCTORS *******************/
            //! Default constructor
            /*
             *  Input arguments
             *  \param time the value to which \c time_ (the internal representation of time) should 
             *  be initialized
             */
            Time (const double& time);
            
            
            /******************* DESTRUCTOR *******************/
            //! Default destructor                
            /*! 
             *  Default destructor, since members of the class are trivially 
             *  destructible.
             */
            virtual ~Time () {};


            /******************* METHODS *******************/
            //! Set time to a given value
            /*!
             *  \param time the value to which \c time_ should be set
             */
            virtual void setTo (const double& time);
            
            //! Get time value
            virtual const double& value ();
            

        // ---------------------------------------------------------------------------------------------//  
        protected:
            //! The time value itself
            double time_;
    };
}

#endif
