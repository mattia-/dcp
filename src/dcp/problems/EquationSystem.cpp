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
    void EquationSystem::solve (const std::string& solveType)
    {
        // this function iterates over solveOrder_ and calls solve (problemName) for each problem, thus delegating
        // to the latter function the task of performing the actual parameters setting and solving

        // check solveType
        if (solveType != "default")
        {
            dolfin::dolfin_error ("dcp: EquationSystem.cpp", 
                                  "solve",
                                  "Unknown solve type \"%s\" requested",
                                  solveType.c_str ());
        }
        dolfin::log (dolfin::DBG, "Selected solve type: %s", solveType.c_str ());

        dolfin::begin (dolfin::DBG, "Start problem solution...");
        
        auto subiterationsBegin = std::find (solveOrder_.begin (), solveOrder_.end (), subiterationsRange_.first);
        auto subiterationsEnd = std::find (solveOrder_.begin (), solveOrder_.end (), subiterationsRange_.second);
        
        auto problemName = solveOrder_.begin ();
        while (problemName != solveOrder_.end ())
        {
            if (problemName != subiterationsBegin)
            {
                dolfin::begin (dolfin::PROGRESS, "Problem: \"%s\"", problemName->c_str ());
                solve_ (*problemName);
                dolfin::end (); // Problem %s

                problemName++;
            }
            else
            {
                subiterate_ (subiterationsBegin, subiterationsEnd);
                problemName = subiterationsEnd;
            }
        }
        
        dolfin::end (); // Solving problems
    }



    /******************* PROTECTED METHODS *******************/
    void EquationSystem::solve_ (const std::string& problemName)
    {
        // get problem with given name from map. 
        dcp::GenericProblem& problem = this -> operator[] (problemName);

        // 1)
        // loop over problemsLinks_ to reset all links to take changes to coefficients 
        // into account. Remember it is a map: elements in it are order according to 
        // the default lexicographic ordering
        dolfin::begin (dolfin::DBG, "Scanning problems links...");

        auto linksIterator = problemsLinks_.begin ();
        while (linksIterator != problemsLinks_.end () && std::get<0> (linksIterator->first) <= problemName)
        {
            if (std::get<0> (linksIterator->first) == problemName)
            {
                linkProblems_ (*linksIterator, storedProblems_);
            }
            ++linksIterator;
        }

        dolfin::end (); // Scanning problems links

        // 2)
        // solve problem
        problem.solve (storedProblemsSolveType_);
    }
}
