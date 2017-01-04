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

#ifndef SRC_PROBLEMS_GENERICLINEARPROBLEM_H_INCLUDE_GUARD
#define SRC_PROBLEMS_GENERICLINEARPROBLEM_H_INCLUDE_GUARD

#include <dcp/problems/GenericProblem.h>
#include <dolfin/la/Matrix.h>
#include <dolfin/la/Vector.h>
#include <dolfin/fem/Form.h>

namespace dcp
{
    /*! \class GenericLinearProblem GenericLinearProblem.h
     *  \brief Class for generic linear differential problems.
     *
     *  This class represents problem of the form
     *  \f[
     *      \mbox{Find } u \in V : a \left(u, v\right) = F \left(v\right) \ \forall\,v\,\in\,V
     *  \f]
     *  with \f$ a \left(u, v\right) : V \times V \rightarrow \mathds{R}\f$ bilinear form on \f$V\f$
     *  and \f$ L \left(v\right) : V \rightarrow \mathds{R} \f$ linear form on the same space.
     *
     *  It inherits publicly from \c GenericProblem  and it extends its functionalities to a generic linear differential
     *  problem.
     *  A concrete linear problem must be declared of type \c dcp::GenericLinearProblem. This class is useful for
     *  polymorphic management of linear problems for which the bilinear and linear forms are not known / accessible in
     *  the current scope.
     */

    class GenericLinearProblem : public dcp::GenericProblem
    {
        // ---------------------------------------------------------------------------------------------//

        public:
            /******************* CONSTRUCTORS *******************/
            //! Default constructor is deleted. The class is not default constructable.
            GenericLinearProblem () = delete;

            //!  Constructor
            /*!
             *  \param functionSpace the problem finite element space
             *  The constructors also sets the following parameters:
             *      - \c "problem_type" a string describing the problem. Default value: \c "linear"
             */
            GenericLinearProblem (const std::shared_ptr<const dolfin::FunctionSpace> functionSpace);


            /******************* DESTRUCTOR *******************/
            //! Destructor
            /*!
             *  Default destructor, since members of the class are trivially
             *  destructible.
             */
            virtual ~GenericLinearProblem () {};


            /******************* GETTERS *******************/
            //! Get const reference to the problem's linear form
            /*!
             *  \return a const reference to the problem's linear form
             */
            virtual const dolfin::Form& bilinearForm () const = 0;

            //! Get const reference to the problem's linear form
            /*!
             *  \return a const reference to the problem's linear form
             */
            virtual const dolfin::Form& linearForm () const = 0;

            //! Get const reference to the problem's linear operator
            /*!
             *  \return a const reference to the problem's linear operator
             */
            virtual const dolfin::GenericMatrix& linearOperator () const = 0;

            //! Get const reference to the problem's right hand side
            /*!
             *  \return a const reference to the problem's right hand side
             */
            virtual const dolfin::GenericVector& rhs () const = 0;


            /******************* METHODS *******************/
            //! Lump system matrix
            /*!
             *  Performs lumping of the system matrix by substituting the diagonal with the sum of each row and
             *  setting the extra-diagonal terms to zero.
             */
            virtual void lumpMatrix () = 0;

            //! Assemble the linear system
            /*!
             *  Assemble the linear system related to the problem to be solved
             */
            virtual void assembleLinearSystem () = 0;

            // ---------------------------------------------------------------------------------------------//

        protected:

            // ---------------------------------------------------------------------------------------------//

        private:
    };
}
#endif
