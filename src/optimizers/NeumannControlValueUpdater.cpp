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

#include <optimizers/NeumannControlValueUpdater.h>
#include <differential_problems/AbstractDifferentialProblem.h>

namespace dcp
{
    /************************* CONSTRUCTORS ********************/
    NeumannControlValueUpdater::NeumannControlValueUpdater (const std::string& problemName, 
                                                            const std::string& coefficientType,
                                                            const std::string& coefficientName) :
        problemName_ (problemName),
        coefficientType_ (coefficientType),
        coefficientName_ (coefficientName)
    {   }



    /************************* OPERATORS ********************/
    void NeumannControlValueUpdater::operator() (dcp::CompositeDifferentialProblem& compositeProblem, 
                                                 const std::shared_ptr <const dolfin::GenericFunction> coefficientValue) const
    {
        dcp::AbstractDifferentialProblem& problem = compositeProblem [problemName_];

        problem.setCoefficient (coefficientType_, coefficientValue, coefficientName_);
    }
}
