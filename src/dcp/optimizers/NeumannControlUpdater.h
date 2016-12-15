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

#ifndef SRC_OPTIMIZERS_NEUMANNCONTROLUPDATER_H_INCLUDE_GUARD
#define SRC_OPTIMIZERS_NEUMANNCONTROLUPDATER_H_INCLUDE_GUARD

#include <dcp/problems/EquationSystem.h>
#include <dolfin/function/GenericFunction.h>
#include <string>
#include <dcp/functions/TimeDependentFunction.h>

namespace dcp
{
    /*! \class NeumannControlUpdater NeumannControlUpdater.h
     *  \brief Class to update the value of the control variable in Neumann boundary control problems.
     *
     *  This class is a functor which can be passed to the constructor of any class of the \c GenericImplementer
     *  hierarchy, which will use it to update the primal problem using the new value of the control parameter as the
     *  optimization proceeds.  Note that only the primal system is updated. If the adjoint system needs to be updated
     *  too a new class should be defined (maybe deriving from this one)
     */
    class NeumannControlUpdater
    {
        // ---------------------------------------------------------------------------------------------//

        public:
            /************************* CONSTRUCTORS ********************/
            //! Default constructor is deleted
            NeumannControlUpdater () = delete;

            //! Constructor
            /*!
             *  Input arguments are:
             *  \param problemName string that identifies the problem (in the \c GenericEquationSystem object
             *  passed as input to <tt>this->operator() ()</tt> ) which contains the control parameter to be updated
             *  \param coefficientType the type of the coefficient representing the control parameter inside the problem.
             *  This will be used by the call to \c dcp::GenericProblem::setCoefficient()
             *  \param coefficientName the name of the coefficient representing the control parameter in the problem
             *  passed as first argument
             */
            NeumannControlUpdater (const std::string& problemName,
                                   const std::string& coefficientType,
                                   const std::string& coefficientName);


            /************************* OPERATORS ********************/
            //! Call operator [1]
            /*!
             *  This will actually perform the updating of the primal problem (which is supposed to be the first element
             *  in \c systems )
             *  Input parametes are:
             *  \param system the system on which to operate
             *  \param coefficientValue the new value for the control parameter identified by \c coefficientName_
             */
            void operator() (const std::vector<std::shared_ptr<dcp::GenericEquationSystem>> systems,
                             const dolfin::GenericFunction& coefficientValue) const;

            //! Call operator [2]
            /*!
             *  This will actually perform the updating of the primal problem (which is supposed to be the first element
             *  in \c systems )
             *  Input parametes are:
             *  \param system the system on which to operate
             *  \param coefficientValue the new value for the control parameter identified by \c coefficientName_
             */
            void operator() (const std::vector<std::shared_ptr<dcp::GenericEquationSystem>> systems,
                             const dcp::TimeDependentFunction& coefficientValue) const;

        // ---------------------------------------------------------------------------------------------//

        protected:
            //! The name identifying the problem inside the \c GenericEquationSystem which contains the
            //! control parameter
            std::string problemName_;

            //! The type of the coefficient representing the control parameter inside the problem. This will be used
            //! by the call to \c dcp::GenericProblem::setCoefficient()
            std::string coefficientType_;

            //! The name of the control parameter inside the problem identified by \c problemName_
            std::string coefficientName_;
    };
}

#endif


