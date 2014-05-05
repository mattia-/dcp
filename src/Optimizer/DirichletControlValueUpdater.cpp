/* 
 *  Copyright (C) 2014, Mattia Tamellini, mattia.tamelllini@gmail.com
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

#include <Optimizer/DirichletControlValueUpdater.hpp>
#include <DifferentialProblem/AbstractDifferentialProblem.hpp>
#include <dolfin/fem/DirichletBC.h>

namespace DCP
{
    /************************* CONSTRUCTORS ********************/
    DirichletControlValueUpdater::DirichletControlValueUpdater (const std::string& problemName, 
                                                                const std::string& dirichletBCName,
                                                                const dolfin::SubDomain& dirichletBoundary,
                                                                boost::shared_ptr<const dolfin::FunctionSpace> functionSpace) : 
        problemName_ (problemName),
        dirichletBCName_ (dirichletBCName),
        dirichletBoundary_ (dirichletBoundary),
        functionSpace_ (functionSpace)
    {  }



    /************************* OPERATORS ********************/
    void DirichletControlValueUpdater::operator() (DCP::CompositeDifferentialProblem& compositeProblem, 
                                                   const dolfin::GenericFunction& dirichletBCValue) const
    {
        DCP::AbstractDifferentialProblem& problem = compositeProblem [problemName_];
        
        problem.removeDirichletBC (dirichletBCName_);
        
        problem.addDirichletBC (dolfin::DirichletBC (*(functionSpace_), dirichletBCValue, dirichletBoundary_),
                                dirichletBCName_);
    }
}
