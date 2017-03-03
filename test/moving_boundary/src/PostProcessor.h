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

#ifndef IVAN_POSTPROCESSOR_H_INCLUDE_GUARD
#define IVAN_POSTPROCESSOR_H_INCLUDE_GUARD

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
#include <dcp/differential_problems/DefaultPostProcessor.h>

#include "balanceTerms.h"

#include "GetPot.h"
extern GetPot inputData;

//TODO: namespace->aegir
namespace Ivan
{

/*! \class PostProcessor PostProcessor.h
 *  \brief Class for output during the solution of dcp::MovingTimeDependentProblem
 *
 *  The problem to be post-processed is stored as a 
 *  <tt> dcp::MovingTimeDependentProblem </tt>
 */

class PostProcessor : public dcp::DefaultPostProcessor
{
    // ---------------------------------------------------------------------------------------------//  

    public:

      PostProcessor (dcp::MovingTimeDependentProblem & pb);
      PostProcessor (dcp::MovingTimeDependentProblem & pb, std::string balanceFileName, std::string intGCLFileName, std::string divGCLFileName, std::string selectedDofsFileName);

      void operator() (int timeStep) { (*this) (timeStep, & (this->pb_)); }
      void operator() (int timeStep, const dcp::MovingTimeDependentProblem * const pb);
      void onOldDomain (int timeStep) { this->onOldDomain (timeStep, & (this->pb_)); }
      void onOldDomain (int timeStep, const dcp::MovingTimeDependentProblem * const pb);

      virtual ~PostProcessor ();

    private:

      std::vector<std::shared_ptr<dolfin::Form> > formsDepOnSol_, formsDepOnOld_, formsDepOnDispl_, formsDepOnDir_;
      balanceTerms::Form_intGCL1::TestSpace gclSpace1_;
      balanceTerms::Form_intGCL2::TestSpace gclSpace2_;
      //std::string outputFileName_, intGCLFileName_, divGCLFileName_, selectedDofsFileName_;
      std::string balanceFileName_, intGCLFileName_, divGCLFileName_, selectedDofsFileName_;
      dolfin::la_index uxDofsNum_, uyDofsNum_, pDofsNum_, wxDofsNum_, wyDofsNum_;
      double * uxDofsVals_, * uyDofsVals_, * pDofsVals_, * wxDofsVals_, * wyDofsVals_;
      const dolfin::la_index * uxDofsIdxs_, * uyDofsIdxs_, * pDofsIdxs_, * wxDofsIdxs_, * wyDofsIdxs_;

//      std::shared_ptr <dcp::AbstractProblem> curvatureProblem_;
      std::shared_ptr<dolfin::Form> curvatureBilinearForm_, curvatureLinearForm_, curvatureAdditionalForm_;
//addNotInMovingAbstract//          const std::vector<dolfin::la_index> & additionalFormDofs_;
//addNotInMovingAbstract//          dolfin::la_index * addIdxs_, * addIdxs_w_;
//addNotInMovingAbstract//          double * addVals_, * addVals_w_;
      dolfin::Matrix curvatureMatrix_;
      dolfin::Vector curvatureTgDiv_, curvatureRhs_;
      dolfin::Function auxiliary_;
      dolfin::LinearSolver linearSolver_;
};

}
#endif
