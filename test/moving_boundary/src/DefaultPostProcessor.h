/* 
 *  Copyright (C) 2017, Ivan Fumagalli, ivan.fumagalli.if@gmail.com
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

#ifndef IVAN_DEFAULTPOSTPROCESSOR_H_INCLUDE_GUARD
#define IVAN_DEFAULTPOSTPROCESSOR_H_INCLUDE_GUARD

#include <vector>
#include <string>
#include <memory>
#include <math.h>

#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshFunction.h>
#include <dolfin/function/Function.h>
#include <dolfin/la/Vector.h>
#include <dolfin/log/dolfin_log.h>
#include <dolfin/parameter/Parameters.h>

#include <dcp/differential_problems/AbstractProblem.h>
#include <dcp/factories/LinearSolverFactory.h>
#include <dcp/differential_problems/SubdomainType.h>
#include <dcp/expressions/TimeDependentExpression.h>
#include <dcp/subdomains/Subdomain.h>
#include <dcp/differential_problems/MeshManager.h>
#include <dcp/differential_problems/utilities.h>

//TODO: namespace->aegir
namespace Ivan
{

class MovingTimeDependentProblem;


/*! \class DefaultPostProcessor DefaultPostProcessor.h
 *  \brief Default interface for output during the solution of MovingTimeDependentProblem
 *
 *  A do-nothing post-processor.
 *  Inherit from this class to pass a post-processor to MovingTimeDependentProblem.
 *
 *  The problem to be post-processed is stored as a 
 *  <tt> Ivan::MovingTimeDependentProblem </tt>
 */

class DefaultPostProcessor
{
    // ---------------------------------------------------------------------------------------------//  

    public:

      DefaultPostProcessor (MovingTimeDependentProblem & pb) :
        pb_ (pb)
      {
      }

      virtual void operator() (int timeStep)
      {
        (*this) (timeStep, & (this->pb_));
      }
      virtual void operator() (int timeStep, const MovingTimeDependentProblem * const pb)
      {
      }

      virtual void onOldDomain (int timeStep)
      {
        this->onOldDomain (timeStep, & (this->pb_));
      }
      virtual void onOldDomain (int timeStep, const MovingTimeDependentProblem * const pb)
      {
      }

      virtual ~DefaultPostProcessor ()
      {
      }

    protected:

      //! The problem whose members should be used
      MovingTimeDependentProblem & pb_;

      //! Coefficients for the forms contained here
      std::map<std::string, std::shared_ptr<dolfin::GenericFunction> > coefficients_;

      //! Forms to be assembled during the post-processing
      std::map<std::string, std::pair<std::shared_ptr<dolfin::Form>, double> > formsNvalues_, formsNvaluesOnOld_;

      //! Forms to be assembled to verify the Geometric Conservation Law
      std::map<std::string, std::shared_ptr<dolfin::Form> > gclForms_;
};

}
#endif
