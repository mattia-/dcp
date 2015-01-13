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

#ifndef SRC_OPTIMIZERS_NEUMANNCONTROLVALUEUPDATER_H_INCLUDE_GUARD
#define SRC_OPTIMIZERS_NEUMANNCONTROLVALUEUPDATER_H_INCLUDE_GUARD

#include <differential_problems/EquationSystem.h>
#include <dolfin/function/GenericFunction.h>
#include <string>

namespace dcp
{
    /*! \class NeumannControlValueUpdater NeumannControlValueUpdater.h
     *  \brief Class to update the value of the control variable in Neumann boundary control problems.
     *  
     *  This class is a functor which can be passed to the method \c apply() of any class
     *  of the \c AbstractOptimizer hierarchy, which will use it to update the value of
     *  the control parameter in the \c EquationSystem (also passed to the method \c apply()
     *  of the same class) as the optimization proceeds.
     */
    class NeumannControlValueUpdater
    {
        // ---------------------------------------------------------------------------------------------//

        public:
            /************************* CONSTRUCTORS ********************/
            //! Default constructor is deleted
            NeumannControlValueUpdater () = delete;
            
            //! Constructor
            /*! 
             *  Input arguments are:
             *  \param problemName string that identifies the problem (in the \c EquationSystem object 
             *  passed as input to <tt>this->operator() ()</tt> ) which contains the control parameter to be updated
             *  \param coefficientType the type of the coefficient representing the control parameter inside the problem. 
             *  This will be used by the call to \c dcp::AbstractProblem::setCoefficient()
             *  \param coefficientName the name of the coefficient representing the control parameter in the problem 
             *  passed as first argument
             */
            NeumannControlValueUpdater (const std::string& problemName, 
                                        const std::string& coefficientType,
                                        const std::string& coefficientName);
            
            
            /************************* OPERATORS ********************/
            //! Call operator. This will actually perform the updating of the control parameter
            /*! 
             *  Input parametes are:
             *  \param compositeProblem the problem on which to operate
             *  \param coefficientValue the new value for the control parameter identified by \c coefficientName_
             */
            void operator() (dcp::EquationSystem& compositeProblem, 
                             const std::shared_ptr <const dolfin::GenericFunction> coefficientValue) const;
 
        // ---------------------------------------------------------------------------------------------//
        
        protected:
            //! The name identifying the problem inside the \c EquationSystem which contains the
            //! control parameter
            std::string problemName_;
           
            //! The type of the coefficient representing the control parameter inside the problem. This will be used
            //! by the call to \c dcp::AbstractProblem::setCoefficient()
            std::string coefficientType_;
           
            //! The name of the control parameter inside the problem identified by \c problemName_
            std::string coefficientName_;
    };
}

#endif

