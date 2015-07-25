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

#include <dcp/problems/EquationSystem.h>
#include <utility>
#include <tuple>
#include <dolfin.h>
#include <algorithm>
#include <dolfin/log/dolfin_log.h>

namespace dcp
{
    /******************* CONSTRUCTORS ******************/
    EquationSystem::EquationSystem () : 
        GenericEquationSystem ()
    { 
        dolfin::log (dolfin::DBG, "EquationSystem object created");
    }

    

    /******************* METHODS *******************/
    void EquationSystem::solve ()
    {
        // this function iterates over solveOrder_ and calls solve (problemName) for each problem, thus delegating
        // to the latter function the task of performing the actual parameters setting and solving
        dolfin::begin ("Solving problems...");
        
        auto subiterationsBegin = std::find (solveOrder_.begin (), solveOrder_.end (), subiterationsRange_.first);
        auto subiterationsEnd = std::find (solveOrder_.begin (), solveOrder_.end (), subiterationsRange_.second);
        
        auto problemName = solveOrder_.begin ();
        while (problemName != solveOrder_.end ())
        {
            if (problemName != subiterationsBegin)
            {
                solve (*problemName);
                problemName++;
            }
            else
            {
                subiterate (subiterationsBegin, subiterationsEnd);
                problemName = subiterationsEnd;
            }
        }
        
        dolfin::end ();
    }



    void EquationSystem::solve (const std::string& problemName)
    {
        dolfin::begin ("Solving problem \"%s\"...", problemName.c_str ());

        // get problem with given name from map. 
        dcp::GenericProblem& problem = this -> operator[] (problemName);

        // 1)
        // loop over problemsLinks_ to reset all links to take changes to coefficients 
        // into account. Remember it is a map: elements in it are order according to 
        // the default lexicographic ordering
        dolfin::begin (dolfin::PROGRESS, "Scanning problems links...");

        auto linksIterator = problemsLinks_.begin ();
        while (linksIterator != problemsLinks_.end () && std::get<0> (linksIterator->first) <= problemName)
        {
            if (std::get<0> (linksIterator->first) == problemName)
            {
                linkProblems (*linksIterator);
            }
            ++linksIterator;
        }

        dolfin::end ();

        // 2)
        // solve problem
        dolfin::log (dolfin::PROGRESS, "Calling solve method on problem...");
        problem.solve (solveType_);
        
        dolfin::end ();
    }
}
