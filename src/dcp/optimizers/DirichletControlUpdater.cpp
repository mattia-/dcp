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

#include <dolfin/fem/DirichletBC.h>
#include <dcp/optimizers/DirichletControlUpdater.h>
#include <dcp/problems/TimeDependentProblem.h>

namespace dcp
{
    /************************* CONSTRUCTORS ********************/
    DirichletControlUpdater::DirichletControlUpdater (const std::string& problemName, 
                                                      const std::string& dirichletBCName,
                                                      const dcp::Subdomain& dirichletBoundary,
                                                      const int& component) : 
        problemName_ (problemName),
        dirichletBCName_ (dirichletBCName),
        dirichletBoundary_ (dirichletBoundary),
        component_ (component)
    {  }



    /************************* OPERATORS ********************/
    void DirichletControlUpdater::operator() (dcp::GenericEquationSystem& system, 
                                              const dolfin::GenericFunction& dirichletBCValue) const
    {
        dcp::GenericProblem& problem = system [problemName_];
        
        problem.removeDirichletBC (dirichletBCName_);
        
        if (component_ < 0)
        {
            problem.addDirichletBC (dolfin::reference_to_no_delete_pointer (dirichletBCValue), 
                                    dolfin::reference_to_no_delete_pointer (dirichletBoundary_),
                                    dirichletBCName_);
        }
        else
        {
            problem.addDirichletBC (dolfin::reference_to_no_delete_pointer (dirichletBCValue), 
                                    dolfin::reference_to_no_delete_pointer (dirichletBoundary_),
                                    component_,
                                    dirichletBCName_);
        }
    }



    void DirichletControlUpdater::operator() (dcp::GenericEquationSystem& system, 
                                              const dcp::TimeDependentFunction& dirichletBCValue) const
    {
        dcp::GenericProblem& problem = system [problemName_];
        auto pointerToProblem = dynamic_cast<dcp::TimeDependentProblem*> (&problem);
        if (pointerToProblem == nullptr)
        {
            dolfin::dolfin_error ("dcp: DirichletControlUpdater.cpp",
                                  "update system",
                                  "Problem \"%s\" in input system is not an object of type dcp::TimeDependentProblem",
                                  problemName_.c_str ());
        }
        
        pointerToProblem->removeTimeDependentDirichletBC (dirichletBCName_);
        
        if (component_ < 0)
        {
            pointerToProblem->addTimeDependentDirichletBC (dirichletBCValue, 
                                                           dirichletBoundary_,
                                                           dirichletBCName_);
        }
        else
        {
            pointerToProblem->addTimeDependentDirichletBC (dolfin::reference_to_no_delete_pointer (dirichletBCValue), 
                                                           dolfin::reference_to_no_delete_pointer (dirichletBoundary_),
                                                           component_,
                                                           dirichletBCName_);
        }
    }
}
