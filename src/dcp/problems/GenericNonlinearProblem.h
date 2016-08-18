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

#ifndef SRC_PROBLEMS_GENERICNONLINEARPROBLEM_H_INCLUDE_GUARD
#define SRC_PROBLEMS_GENERICNONLINEARPROBLEM_H_INCLUDE_GUARD

#include <dcp/problems/GenericProblem.h>
#include <dolfin/fem/Form.h>

namespace dcp
{
    /*! \class GenericNonlinearProblem NonlinearProblem.h
     *  \brief Generic class for non-linear differential problems.
     *
     *  This class represents problem of the form
     *  \f[
     *      \mbox{Find } u \in V : F \left(u, v\right) = 0 \ \forall\,v\,\in\,V
     *  \f]
     *  with \f$ F \left(u, v\right) : V \times V \rightarrow \mathds{R}\f$ non linear in the problem unknown \f$u\f$.
     *  
     *  It inherits publicly from \c GenericProblem and it extends its functionalities to a concrete differential
     *  problem.
     *  A concrete nonlinear problem must be declared of type \c dcp::GenericLinearProblem. This class is useful for 
     *  polymorphic management of nonlinear problems for which the residual and jacobian forms are not known / accessible 
     *  in the current scope.
     */

    class GenericNonlinearProblem : public dcp::GenericProblem
    {
        // ---------------------------------------------------------------------------------------------//  

        public:
            /******************* CONSTRUCTORS *******************/
            //! Default constructor is deleted. The class is not default constructable.
            GenericNonlinearProblem () = delete;

            //!  Constructor
            /*!
             *  \param functionSpace the problem finite element space
             *  The constructors also sets the following parameters:
             *      - \c "problem_type" a string describing the problem. Default value: \c "non_linear"
             */
            GenericNonlinearProblem (const std::shared_ptr<const dolfin::FunctionSpace> functionSpace);


            /******************* DESTRUCTOR *******************/
            //! Destructor
            /*! Default destructor, since members of the class are trivially 
             * destructible.
             */
            virtual ~GenericNonlinearProblem () {};


            /******************* GETTERS *******************/
            //! Get const reference to the problem's residual form
            /*! 
             *  \return a const reference to the problem's residual form
             */
            virtual const dolfin::Form& residualForm () const = 0;

            //! Get const reference to the problem's jacobian form
            /*! 
             *  \return a const reference to the problem's jacobian form
             */
            virtual const dolfin::Form& jacobianForm () const = 0;

            // ---------------------------------------------------------------------------------------------//

        protected:

            // ---------------------------------------------------------------------------------------------//

        private:
    };
}
#endif

