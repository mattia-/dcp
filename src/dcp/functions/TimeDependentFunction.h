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

#ifndef SRC_FUNCTIONS_TIMEDEPENDENTFUNCTION_H_INCLUDE_GUARD
#define SRC_FUNCTIONS_TIMEDEPENDENTFUNCTION_H_INCLUDE_GUARD

#include <dolfin/function/FunctionSpace.h>
#include <dolfin/function/Function.h>
#include <dolfin/io/File.h>
#include <dolfin/plot/VTKPlotter.h>
#include <vector>
#include <utility>

namespace dcp
{
    /*! \class TimeDependentFunction TimeDependentFunction.h
     *  \brief Class holding the solution of problems in the \c dcp hierarchy
     *
     *  This class contains the interface for the solution of problems in the \c dcp hierarchy. It derives from
     *  \c std::vector<std::pair<double,dolfin::Function>> and extends its functionalities to some needed special cases.
     *
     *  IMPORTANT: even though no control is in place, times in dcp::TimeDependentFunction objects are always supposed
     *  to be sorted in ascending order within the container
     */
    class TimeDependentFunction : public std::vector<std::pair<double, dolfin::Function>>
    {
        // ---------------------------------------------------------------------------------------------//

        public:
            /************************* CONSTRUCTORS ********************/
            //! Default constructor
            TimeDependentFunction ();

            //! Fill constructor [1]
            /*!
             *  Creates object with given size containing pairs made from the input initial time \c t0 (incremented by
             *  \c dt for each element) and functions build on the input function space.
             *
             *  \param n the size
             *  \param t0 the initial time
             *  \param dt the time step
             *  \param functionSpace the function space
             */
            TimeDependentFunction (const int &n,
                                   const double& t0,
                                   const double& dt,
                                   const std::shared_ptr<const dolfin::FunctionSpace> functionSpace);

            //! Fill constructor [2]
            /*!
             *  Creates object containing pairs made from the input initial time \c t0 incremented by \c dt for each
             *  element until final time is reached and functions build on the input function space.
             *
             *  \param t0 the initial time
             *  \param dt the time step
             *  \param tf the final time
             *  \param functionSpace the function space
             */
            TimeDependentFunction (const double& t0,
                                   const double& dt,
                                   const double& tf,
                                   const std::shared_ptr<const dolfin::FunctionSpace> functionSpace);


            //! Fill constructor [3]
            /*!
             *  Creates object with given size containing pairs made from the input initial time \c t0 (incremented by
             *  \c dt for each element) and copies of the input funciton
             *
             *  \param n the size
             *  \param t0 the initial time
             *  \param dt the time step
             *  \param function the function to be copied in the container
             */
            TimeDependentFunction (const int &n,
                                   const double& t0,
                                   const double& dt,
                                   const dolfin::Function& function);

            //! Fill constructor [2]
            /*!
             *  Creates object containing pairs made from the input initial time \c t0 incremented by \c dt for each
             *  element until final time is reached and copies of the input function
             *
             *  \param t0 the initial time
             *  \param dt the time step
             *  \param tf the final time
             *  \param function the function to be copied in the container
             */
            TimeDependentFunction (const double& t0,
                                   const double& dt,
                                   const double& tf,
                                   const dolfin::Function& function);

            //! Fill constructor [5]
            /*!
             *  Creates object with given size containing copies of the input pair
             *  input function space.
             *
             *  \param n the size
             *  \param value the pair to be stored in the container
             */
            TimeDependentFunction (const int &n,
                                   const TimeDependentFunction::value_type& value);

            //! Copy constructor
            /*!
             *  \param other the time-dependent function to copy to build the new object
             */
            TimeDependentFunction (const TimeDependentFunction& other) = default;


            /************************* DESTRUCTOR ********************/
            //! Destructor
            virtual ~TimeDependentFunction () {};


            /********************** METHODS ***********************/
            //! Plot time dependent function
            /*!
             *  \param title the title of the plot; time value will be appended to it
             *  \param pause switch to enable (\c true) or disable (\c false) pausing after each plot
             */
            void plot (std::string title, const bool& pause);

            //! Compute time integral of the function
            /*!
             *  The time integral is a function in the same function space which at each point has value equal to the
             *  integral in time of the function at that same point. If the time dependent function contains no data, an
             *  error is issued. If instead the time dependent function only has one element, zero (as a function) is
             *  returned.
             *
             *  \return the result of the time integration (move semantic is used if implemented by FEniCS)
             */
            dolfin::Function computeTimeIntegral () const;

            /********************** OPERATORS ***********************/
            //! Assignement operator
            /*!
             *  This allows to assign values to a time-dependent function from an expression, by assigning the
             *  evaluated expression at each timestep to each element in the time-dependent function.
             *  If the expression is constant with respect to time, the time-dependent function will simply contain
             *  the same values at each timestep. If the expression is indeed time-dependent, at the end of the function
             *  the value of time contained in the expression at the function call will be restored.
             *
             *  \param expression the expression;
             */
            dcp::TimeDependentFunction& operator= (dolfin::Expression& expression);

            //! Sum operator
            /*!
             *  Performs the sum <tt>a + b</tt>.
             *  The sum is only performed if the time values in the two functions match. In this case, the normal sum
             *  between \c dolfin::Function is performed at each timestep. An error is issued otherwise.
             *  \param left the first function of the sum
             *  \param right the second function of the sum
             *
             *  Reminder to developers: C++11 already implements move semantic or return-value optimization if possible
             *  when returning an object from a function, so there is no need to explicitly use the move semantic
             */
            friend dcp::TimeDependentFunction operator+ (const dcp::TimeDependentFunction& left,
                                                         const dcp::TimeDependentFunction& right);

            //! Difference operator
            /*!
             *  Performs the difference <tt>a - b</tt>.
             *  The sum is only performed if the time values in the two functions match. In this case, the normal
             *  difference between \c dolfin::Function is performed at each timestep. An error is issued otherwise.
             *  \param left the first function of the difference
             *  \param right the second function of the difference
             *
             *  Reminder to developers: C++11 already implements move semantic or return-value optimization if possible
             *  when returning an object from a function, so there is no need to explicitly use the move semantic
             */
            friend dcp::TimeDependentFunction operator- (const dcp::TimeDependentFunction& left,
                                                         const dcp::TimeDependentFunction& right);

            //! Multiplication by scalar operator [1]
            /*!
             *  Performs the multiplication <tt>a * k</tt> where \c a is of type \c dcp::TimeDependentFunction and \c k
             *  is a \c double .
             *  The normal multiplication between \c dolfin::Function and \c double is performed at each timestep.
             *  \param function the function
             *  \param scalar the scalar
             *
             *  Reminder to developers: C++11 already implements move semantic or return-value optimization if possible
             *  when returning an object from a function, so there is no need to explicitly use the move semantic
             */
            friend dcp::TimeDependentFunction operator* (const dcp::TimeDependentFunction& function,
                                                         const double& scalar);

            //! Multiplication by scalar operator [2]
            /*!
             *  Performs the multiplication <tt>k * a</tt> where \c k is a \c double and \c a is of type
             *  \c dcp::TimeDependentFunction .
             *  The normal multiplication between \c dolfin::Function and \c double is performed at each timestep.
             *  \param scalar the scalar
             *  \param function the function
             *
             *  Reminder to developers: C++11 already implements move semantic or return-value optimization if possible
             *  when returning an object from a function, so there is no need to explicitly use the move semantic
             */
            friend dcp::TimeDependentFunction operator* (const double& scalar,
                                                         const dcp::TimeDependentFunction& function);

            //! Multiplication by scalar operator [3]
            /*!
             *  Performs the multiplication <tt>a * k</tt> where \c a is of type \c dcp::TimeDependentFunction and \c k
             *  is a \c int .
             *  The normal multiplication between \c dolfin::Function and \c int is performed at each timestep.
             *  \param function the function
             *  \param scalar the scalar
             *
             *  Reminder to developers: C++11 already implements move semantic or return-value optimization if possible
             *  when returning an object from a function, so there is no need to explicitly use the move semantic
             */
            friend dcp::TimeDependentFunction operator* (const dcp::TimeDependentFunction& function,
                                                         const int& scalar);

            //! Multiplication by scalar operator [4]
            /*!
             *  Performs the multiplication <tt>k * a</tt> where \c k is a \c int and \c a is of type
             *  \c dcp::TimeDependentFunction .
             *  The normal multiplication between \c dolfin::Function and \c int is performed at each timestep.
             *  \param scalar the scalar
             *  \param function the function
             *
             *  Reminder to developers: C++11 already implements move semantic or return-value optimization if possible
             *  when returning an object from a function, so there is no need to explicitly use the move semantic
             */
            friend dcp::TimeDependentFunction operator* (const int& scalar,
                                                         const dcp::TimeDependentFunction& function);

            //! Unary minus operator
            /*!
             *  Computes the function defined as <tt>-a</tt>. It simply multiplies the given time-dependent function
             *  by \c -1 .
             *  \param function the function
             *
             *  Reminder to developers: C++11 already implements move semantic or return-value optimization if possible
             *  when returning an object from a function, so there is no need to explicitly use the move semantic
             */
            friend dcp::TimeDependentFunction operator- (const dcp::TimeDependentFunction& function);

            //! Write to file operator
            /*!
             *  Writes to the given \c dolfin::file the given time dependent function
             *
             *  \param file the file
             *  \param function the function
             */
            friend void operator<< (dolfin::File& file, const dcp::TimeDependentFunction& function);

            /********************** VARIABLES ***********************/
            //! the function parameters; may be useful in derived classes
            dolfin::Parameters parameters;

            // ---------------------------------------------------------------------------------------------//

        protected:
            //! The plotter for the function
            std::shared_ptr<dolfin::VTKPlotter> plotter_;


            // ---------------------------------------------------------------------------------------------//

        private:

    };

}
#endif

