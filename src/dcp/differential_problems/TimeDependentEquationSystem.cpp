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
    void TimeDependentEquationSystem::addLinkToPreviousSolution (const std::string& linkFrom, 
                                                                 const std::string& linkedCoefficientName,
                                                                 const std::string& linkedCoefficientType, 
                                                                 const std::string& linkTo,
                                                                 const int& nStepsBack,
                                                                 const bool& forceRelinking)
    {
        dolfin::begin (dolfin::DBG, 
                       "Setting up link (%s, %s, %s) -> (%s, all solution components, %d time steps back)...",
                       linkFrom.c_str (),
                       linkedCoefficientName.c_str (),
                       linkedCoefficientType.c_str (),
                       linkTo.c_str (),
                       nStepsBack);
        
        // create pair containing the link information passed as input arguments.
        auto link = std::make_pair (std::make_tuple (linkFrom, linkedCoefficientName, linkedCoefficientType), 
                                    std::make_tuple (linkTo, -1, nStepsBack));

        // search for map key in linksToPreviousSolutions_. 
        auto linkPosition = linksToPreviousSolutions_.find (link.first);

        if (linkPosition == linksToPreviousSolutions_.end ()) // if key not found in map, insert link
        {
            dolfin::log (dolfin::DBG, "Inserting link in links to previous solutions map...");
            linksToPreviousSolutions_.insert (link);
            
            // perform linking
            linkProblemToPreviousSolution (link);
        }
        else if (forceRelinking == true) // if key found in map but forceRelinking set to true, erase 
        // current link and insert the new one
        {
            dolfin::cout << "In equation system: erasing link:" << dolfin::endl;
            dolfin::cout << "\t(" 
                << std::get<0> (linkPosition->first) 
                << ", " 
                << std::get<1> (linkPosition->first) 
                << ", " 
                << std::get<2> (linkPosition->first) 
                << ") -> (" 
                << std::get<0> (linkPosition->second)
                << ", "
                << std::string (std::get<1> (linkPosition->second) == -1 ? 
                                "all solution components, " : 
                                "component " + std::to_string (std::get<1> (linkPosition->second)) + ", ")
                << std::get<2> (linkPosition->second)
                << " time steps back)"
                << dolfin::endl;

            linksToPreviousSolutions_.erase (linkPosition);

            dolfin::cout << "and inserting link: " << dolfin::endl;
            dolfin::cout << "\t(" 
                << std::get<0> (link.first)
                << ", " 
                << std::get<1> (link.first)
                << ", " 
                << std::get<2> (link.first)
                << ") -> (" 
                << std::get<0> (link.second)
                << ", all solution components, " 
                << std::get<2> (link.second)
                << " time steps back)"
                << dolfin::endl;

            linksToPreviousSolutions_.insert (link);
            
            // perform linking
            linkProblemToPreviousSolution (link);
        }
        else
        {
            dolfin::warning 
                ("link (%s, %s, %s) -> (%s, all solution components, %d time steps back) not added. Key is already present in map",
                 linkFrom.c_str (),
                 linkedCoefficientName.c_str (),
                 linkedCoefficientType.c_str (),
                 linkTo.c_str (),
                 nStepsBack);
        }
        dolfin::end ();
    }



    void TimeDependentEquationSystem::addLinkToPreviousSolution (const std::string& linkFrom, 
                                                                 const std::string& linkedCoefficientName,
                                                                 const std::string& linkedCoefficientType, 
                                                                 const std::string& linkTo,
                                                                 const int& linkToComponent,
                                                                 const int& nStepsBack,
                                                                 const bool& forceRelinking) 
    {
        dolfin::begin (dolfin::DBG, 
                       "Setting up link (%s, %s, %s) -> (%s, component %d, %d time steps back)...",
                       linkFrom.c_str (),
                       linkedCoefficientName.c_str (),
                       linkedCoefficientType.c_str (),
                       linkTo.c_str (),
                       linkToComponent,
                       nStepsBack);
        
        // create pair containing the link information passed as input arguments.
        auto link = std::make_pair (std::make_tuple (linkFrom, linkedCoefficientName, linkedCoefficientType), 
                                    std::make_tuple (linkTo, linkToComponent, nStepsBack));

        // search for map key in linksToPreviousSolutions_. 
        auto linkPosition = linksToPreviousSolutions_.find (link.first);

        if (linkPosition == linksToPreviousSolutions_.end ()) // if key not found in map, insert link
        {
            dolfin::log (dolfin::DBG, "Inserting link in links map...");
            linksToPreviousSolutions_.insert (link);
            
            // perform linking
            linkProblemToPreviousSolution (link);
        }
        else if (forceRelinking == true) // if key found in map but forceRelinking set to true, erase 
        // current link and insert the new one
        {
            dolfin::cout << "In equation system: erasing link:" << dolfin::endl;
            dolfin::cout << "\t(" 
                << std::get<0> (linkPosition->first) 
                << ", " 
                << std::get<1> (linkPosition->first) 
                << ", " 
                << std::get<2> (linkPosition->first) 
                << ") -> (" 
                << std::get<0> (linkPosition->second)
                << ", "
                << std::string (std::get<1> (linkPosition->second) == -1 ? 
                                "all solution components, " : 
                                "component " + std::to_string (std::get<1> (linkPosition->second)) + ", ")
                << std::get<2> (linkPosition->second)
                << " time steps back)"
                << dolfin::endl;

            linksToPreviousSolutions_.erase (linkPosition);

            dolfin::cout << "and inserting link: " << dolfin::endl;
            dolfin::cout << "\t(" 
                << std::get<0> (link.first)
                << ", " 
                << std::get<1> (link.first)
                << ", " 
                << std::get<2> (link.first)
                << ") -> (" 
                << std::get<0> (link.second)
                << ", component " 
                << std::get<1> (link.second)
                << ", "
                << std::get<2> (link.second)
                << " time steps back)"
                << dolfin::endl;

            linksToPreviousSolutions_.insert (link);
            
            // perform linking
            linkProblemToPreviousSolution (link);
        }
        else
        {
            dolfin::warning 
                ("link (%s, %s, %s) -> (%s, component %d, %d time steps back) not added. Key is already present in map",
                 linkFrom.c_str (),
                 linkedCoefficientName.c_str (),
                 linkedCoefficientType.c_str (),
                 linkTo.c_str (),
                 linkToComponent,
                 nStepsBack);
        }
        dolfin::end ();
    }
            


    bool TimeDependentEquationSystem::removeLinkToPreviousSolution (const std::string& linkFrom, 
                                                                    const std::string& linkedCoefficientName,
                                                                    const std::string& linkedCoefficientType)
    {
        auto linkKey = std::make_tuple (linkFrom, linkedCoefficientName, linkedCoefficientType);
        return (linksToPreviousSolutions_.erase (linkKey)) == 1 ? true : false;
    }
            


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
    


    const dcp::TimeDependentProblem& TimeDependentEquationSystem::operator[] (const std::string& name) const
    {
        auto problemIterator = storedProblems_.find (name);
        if (problemIterator == storedProblems_.end ())
        {
            dolfin::dolfin_error ("dcp: TimeDependentEquationSystem.cpp",
                                  "operator[]", 
                                  "Problem \"%s\" not found in stored problems map", 
                                  name.c_str ());
        }
        return *(std::dynamic_pointer_cast<dcp::TimeDependentProblem> (problemIterator->second));
    }



    dcp::TimeDependentProblem& TimeDependentEquationSystem::operator[] (const std::string& name)
    {
        auto problemIterator = storedProblems_.find (name);
        if (problemIterator == storedProblems_.end ())
        {
            dolfin::dolfin_error ("dcp: TimeDependentEquationSystem.cpp",
                                  "operator[]", 
                                  "Problem \"%s\" not found in stored problems map", 
                                  name.c_str ());
        }
        return *(std::dynamic_pointer_cast<dcp::TimeDependentProblem> (problemIterator->second));
    }



    const dcp::TimeDependentProblem& TimeDependentEquationSystem::operator[] (const std::size_t& position) const
    {
        if (position >= solveOrder_.size ())
        {
            dolfin::dolfin_error ("dcp: TimeDependentEquationSystem.cpp",
                                  "operator[]",
                                  "Input value \"%d\" is greater than problems vector size",
                                  position);
        }
        return this->operator[] (solveOrder_ [position]);
    }



    dcp::TimeDependentProblem& TimeDependentEquationSystem::operator[] (const std::size_t& position)
    {
        if (position >= solveOrder_.size ())
        {
            dolfin::dolfin_error ("dcp: TimeDependentEquationSystem.cpp",
                                  "operator[]",
                                  "Input value \"%d\" is greater than problems vector size",
                                  position);
        }
        return this->operator[] (solveOrder_ [position]);
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
            
            auto previousSolutionsLinksIterator = linksToPreviousSolutions_.begin ();
            while (previousSolutionsLinksIterator != linksToPreviousSolutions_.end () 
                   && 
                   std::get<0> (previousSolutionsLinksIterator->first) <= problemName)
            {
                if (std::get<0> (previousSolutionsLinksIterator->first) == problemName)
                {
                    linkProblemToPreviousSolution (*previousSolutionsLinksIterator);
                }
                ++previousSolutionsLinksIterator;
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
    


    /******************* PROTECTED METHODS *******************/
    void TimeDependentEquationSystem::linkProblemToPreviousSolution (const PreviousSolutionLink& link)
    {
        if (std::get<1> (link.second) == -1)
        {
            dolfin::begin (dolfin::DBG, 
                           "Considering link: (%s, %s, %s) -> (%s, all solution componentes, %d time steps back)...",
                           (std::get<0> (link.first)).c_str (),
                           (std::get<1> (link.first)).c_str (),
                           (std::get<2> (link.first)).c_str (),
                           (std::get<0> (link.second)).c_str (),
                            std::get<2> (link.second));
        }
        else
        {
            dolfin::begin (dolfin::DBG, 
                           "Considering link: (%s, %s, %s) -> (%s, component %d, %d time steps back)...",
                           (std::get<0> (link.first)).c_str (),
                           (std::get<1> (link.first)).c_str (),
                           (std::get<2> (link.first)).c_str (),
                           (std::get<0> (link.second)).c_str (),
                            std::get<1> (link.second),
                            std::get<2> (link.second));
        }
        
        // check if problem that needs linking exists
        dolfin::log (dolfin::DBG, 
                     "Looking for problem \"%s\" in problems map...", 
                     (std::get<0> (link.first)).c_str ());
        auto problemIterator = storedProblems_.find (std::get<0> (link.first));

        if (problemIterator == storedProblems_.end ())
        {
            dolfin::warning ("Problem \"%s\" not found in stored problems map", 
                             (std::get<0> (link.first)).c_str ());
            dolfin::end ();
            return;
        }

        dcp::AbstractProblem& problem = *(problemIterator->second);

        // check if target problem of the link exists
        dolfin::log (dolfin::DBG, "Looking for link target in problems map...");
        auto targetProblemIterator = storedProblems_.find (std::get<0> (link.second));
        if (targetProblemIterator == storedProblems_.end ())
        {
            dolfin::warning ("Cannot link problem \"%s\". No such problem found in stored problems map",
                             (std::get<0> (link.first)).c_str ());
            dolfin::end ();
            return;
        }

        // unlike in AbstractEquationSystem, we use a reference to TimeDependentProblem (i.e. the derived class).
        // Two reasons for this:
        // 1) we need to call solutionsVector
        // 2) if it is not a TimeDependentProblem, the whole method does not make sense
        dcp::TimeDependentProblem& targetProblem = 
            static_cast<dcp::TimeDependentProblem&> (*(targetProblemIterator->second));

        if (std::get<1> (link.second) == -1)
        {
            dolfin::log 
                (dolfin::DBG, 
                 "Linking coefficient \"%s\" of type \"%s\" of problem \"%s\" to solution of problem \"%s\" from %d time steps back...",
                 (std::get<1> (link.first)).c_str (),
                 (std::get<2> (link.first)).c_str (),
                 (std::get<0> (link.first)).c_str (),
                 (std::get<0> (link.second)).c_str (),
                  std::get<2> (link.second));
            
            // get target problem solution vector
            auto& targetProblemSolutionsVector = targetProblem.solutions ();  
            
            // get target function, by going back from the last element of nStepsBack steps. 
            // NB: we use operator+ to traverse the vector backwards, since rbegin is a REVERSE iterator
            int nStepsBack = std::get<2> (link.second);
            const dolfin::Function& targetFunction = (targetProblemSolutionsVector.rbegin() + nStepsBack)->second;

            problem.setCoefficient (std::get<2> (link.first), 
                                    dolfin::reference_to_no_delete_pointer (targetFunction),
                                    std::get<1> (link.first));
        }
        else
        {
            dolfin::log 
                (dolfin::DBG, 
                 "Linking coefficient \"%s\" of type \"%s\" of problem \"%s\" to component %d of solution of problem \"%s\" from %d time steps back...",
                 (std::get<1> (link.first)).c_str (),
                 (std::get<2> (link.first)).c_str (),
                 (std::get<0> (link.first)).c_str (),
                 std::get<1> (link.second),
                 (std::get<0> (link.second)).c_str (),
                 std::get<2> (link.second));

            // get target problem solution vector
            auto& targetProblemSolutionsVector = targetProblem.solutions ();  
            
            // get target function, by going back from the last element of nStepsBack steps. 
            // NB: we use operator+ to traverse the vector backwards, since rbegin is a REVERSE iterator
            int nStepsBack = std::get<2> (link.second);
            const dolfin::Function& targetFunction = (targetProblemSolutionsVector.rbegin() + nStepsBack)->second;

            int component = std::get<1> (link.second);
            problem.setCoefficient (std::get<2> (link.first), 
                                    dolfin::reference_to_no_delete_pointer (targetFunction [component]),
                                    std::get<1> (link.first));

        }
        
        dolfin::end ();
    }
}
