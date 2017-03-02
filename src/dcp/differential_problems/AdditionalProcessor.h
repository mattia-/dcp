/* 
 *  Copyright (C) 2015, Ivan Fumagalli, ivan.fumagalli@polimi.it
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

#ifndef DCP_DIFFERENTIAL_PROBLEMS_ADDITIONALPROCESSOR_H_INCLUDE_GUARD
#define DCP_DIFFERENTIAL_PROBLEMS_ADDITIONALPROCESSOR_H_INCLUDE_GUARD

#include <dolfin.h>

namespace dcp
{

/*! \class AdditionalProcessor AdditionalProcessor.h
 *  \brief Functor giving the general interface for the processing of additional terms to be given to #MovingLinearProblem
 *
 * Inherit from this class to define a processor to pass to the ctor of #MovingLinearProblem.
 * Otherwise, it is a do-nothing functor.
 */
class AdditionalProcessor
{
  public:

    virtual void operator() (dolfin::Vector & processedVec, const dolfin::Vector & vec, const std::vector<dolfin::la_index> & additionalDofs)
    {
    }
};

} //end of namespace
#endif
