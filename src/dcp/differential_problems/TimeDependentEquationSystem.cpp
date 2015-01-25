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

#include <dcp/differential_problems/TimeDependentEquationSystem.h>
#include <dcp/differential_problems/TimeDependentProblem.h>
#include <utility>
#include <tuple>
#include <dolfin.h>
#include <algorithm>
#include <dolfin/log/dolfin_log.h>

namespace dcp
{
    /******************* CONSTRUCTORS ******************/
    TimeDependentEquationSystem::TimeDependentEquationSystem () : 
        AbstractEquationSystem ()
    { 
        dolfin::log (dolfin::DBG, "TimeDependentEquationSystem object created");
    }

    

    /******************* METHODS *******************/
    bool TimeDependentEquationSystem::isFinished ()
    {
        std::size_t nFinished = 0;
        for (auto mapElement : storedProblems_)
        {
            dcp::TimeDependentProblem& tmp = static_cast<dcp::TimeDependentProblem&> (*(mapElement.second));
            nFinished += tmp.isFinished ();
        }
        return (nFinished == storedProblems_.size ());
    }
    


    void TimeDependentEquationSystem::solve (const bool& forceRelinking)
    {
        // this function iterates over solveOrder_ and calls solve (problemName) for each problem, thus delegating
        // to the latter function the task of performing the actual parameters setting and solving.
        // The loop is repeated until isFinished() returns true, that is until all problems' time loops are ended
        dolfin::begin ("Solving problems...");
        
        int iterationCounter = 0;
        while (isFinished () == 0)
        {
            iterationCounter++;
            dolfin::log (dolfin::INFO, "=====================================");
            dolfin::log (dolfin::INFO, "TIME DEPENDENT SYSTEM ITERATION: %d", iterationCounter);
            dolfin::begin (dolfin::INFO, "=====================================");
            for (auto problem : solveOrder_)
            {
                solve (problem, forceRelinking);
            }
            dolfin::log (dolfin::INFO, "");
            dolfin::end ();
        }
        
        dolfin::end ();
    }



    void TimeDependentEquationSystem::solve (const std::string& problemName, const bool& forceRelinking)
    {
        dolfin::begin ("Solving problem \"%s\"...", problemName.c_str ());

        // get problem with given name from map. Variable problemIterator will be a
        // std::map <std::string, std::unique_ptr <dcp::AbstractProblem>::iterator
        dolfin::log (dolfin::DBG, "Looking for problem \"%s\" in problems map...", problemName.c_str ());
        auto problemIterator = storedProblems_.find (problemName);

        if (problemIterator == storedProblems_.end ())
        {
            dolfin::warning ("Problem \"%s\" not found in stored problems map", problemName.c_str ());
            return;
        }

        dcp::AbstractProblem& problem = *(problemIterator->second);

        // 1)
        // if forceRelinking is true, loop over problemsLinks_. 
        // Remember it is a map. Elements in it are order according to the default
        // lexicographical ordering
        if (forceRelinking == true)
        {
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
        }

        // 2)
        // solve problem
        dolfin::log (dolfin::PROGRESS, "Calling solve method on problem...");
        problem.solve ("step");
        
        dolfin::end ();
    }



    void TimeDependentEquationSystem::solve (const char* problemName, const bool& forceRelinking)
    {
        solve (std::string (problemName), forceRelinking);
    }
}
