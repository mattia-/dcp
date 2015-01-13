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

#include <differential_problems/AbstractEquationSystem.h>
#include <utility>
#include <tuple>
#include <dolfin.h>
#include <algorithm>
#include <dolfin/log/dolfin_log.h>

namespace dcp
{
    /******************* CONSTRUCTORS ******************/
    AbstractEquationSystem::AbstractEquationSystem () : 
        storedProblems_ (),
        solveOrder_ (),
        problemsLinks_ ()
    { 
        dolfin::log (dolfin::DBG, "AbstractEquationSystem object created");
    }

    

    /******************* METHODS *******************/
    std::size_t AbstractEquationSystem::size ()
    {
        return storedProblems_.size ();
    }



    void AbstractEquationSystem::addProblem (const std::string& problemName, 
                                             dcp::AbstractProblem& problem)
    {
        dolfin::begin (dolfin::DBG, 
                       "Inserting problem \"%s\" in equation system...", 
                       problemName.c_str ());
         
        // insert problem into storedProblems_ 
        dolfin::log (dolfin::DBG, 
                     "Inserting problem in problems map with name \"%s\"...", 
                     problemName.c_str ());
        auto result = storedProblems_.insert 
            (std::make_pair (problemName, dolfin::reference_to_no_delete_pointer (problem)));
        
        // if problem was not inserted in list, issue a warning; else, add it also to the vector solveOrder_
        if (result.second == false)
        {
            dolfin::warning ("Problem \"%s\" already exist in equation system", 
                             problemName.c_str ());
        }
        else
        {
            dolfin::log (dolfin::DBG, 
                         "Inserting problem in problem-names vector with name \"%s\"...", 
                         problemName.c_str ());
            solveOrder_.emplace_back (problemName);
        }
        dolfin::end ();
    }
    
    

    void AbstractEquationSystem::addProblem (const std::string& problemName, 
                                             const std::shared_ptr<dcp::AbstractProblem> problem)
    {
        dolfin::begin (dolfin::DBG, 
                       "Inserting problem \"%s\" in equation system...", 
                       problemName.c_str ());
        
        // insert problem into storedProblems_ 
        dolfin::log (dolfin::DBG, 
                     "Inserting problem in problems map with name \"%s\"...", 
                     problemName.c_str ());
        auto result = storedProblems_.insert (std::make_pair (problemName, problem));
        
        // if problem was not inserted in list, issue a warning; else, add it also to the vector solveOrder_
        if (result.second == false)
        {
            dolfin::warning ("Problem \"%s\" already exist in equation system", 
                             problemName.c_str ());
        }
        else
        {
            dolfin::log (dolfin::DBG, 
                         "Inserting problem in problem-names vector with name \"%s\"...", 
                         problemName.c_str ());
            solveOrder_.emplace_back (problemName);
        }
        dolfin::end ();
    }



    void AbstractEquationSystem::removeProblem (const std::string& problemName)
    {
        dolfin::begin (dolfin::DBG, 
                       "Removing problem \"%s\" from equation system...", 
                       problemName.c_str ());
        
        // delete problem from storedProblems_
        dolfin::log (dolfin::DBG, 
                     "Removing problem \"%s\" from problems map...", 
                     problemName.c_str ());
        auto result = storedProblems_.erase (problemName);
        
        // if problem was not inserted in list, issue a warning; else, erase it also from solveOrder_
        // remember that erase returns the number of elements removed, which in the case of a map is at most 1
        if (result < 1) 
        {
            dolfin::warning ("Problem \"%s\" was not found in equation system. Maybe you used a wrong name?",
                             problemName.c_str ());
        }
        else
        {
            // 1)
            // remove problemName from solveOrder_. Remember that we cannot be sure that such name is 
            // unique in vector and that problems are not in any kind of order.
            // std::find will return an iterator to the element if found and vector::end () if not found
            dolfin::log (dolfin::DBG, 
                         "Removing every occurrence of problem \"%s\" from problem-names vector...", 
                         problemName.c_str ());
            int erasedCount = 0;
            auto problemPosition = find (solveOrder_.begin (), solveOrder_.end (), problemName);
            while (problemPosition != solveOrder_.end ())
            {
                ++erasedCount;
                solveOrder_.erase (problemPosition);
                problemPosition = find (solveOrder_.begin (), solveOrder_.end (), problemName);
            }
            dolfin::log (dolfin::DBG, 
                         "Removed %d entries from problem-names vector", 
                         erasedCount);
            
            // 2)
            // remove problemName from problemsLinks_.
            // We use an important statement from the c++ standard: 
            // When erasing from a map, iterators, pointers and references referring to elements removed by the 
            // function are invalidated. All other iterators, pointers and references keep their validity.
            dolfin::log (dolfin::DBG, 
                         "Removing every occurrence of problem \"%s\" from links map...", 
                         problemName.c_str ());
            erasedCount = 0;
            auto linksIterator = problemsLinks_.begin (); // iterator pointing to the first element of the set
            while (linksIterator != problemsLinks_.end ())
            {
                // delete element if problemName appears either as first or as fourth string in the map
                if (std::get<0> (linksIterator->first) == problemName 
                    || 
                    std::get<0> (linksIterator->second) == problemName)
                {
                    auto auxIterator = linksIterator; // this will be used for the call to function erase
                    ++linksIterator;   // iterator incremented to point to next element. 
                                       // This is performed before erasing element, 
                                       // so that increment is still valid
                    
                    problemsLinks_.erase (auxIterator);
                    ++erasedCount;
                }
                else
                {
                    ++linksIterator;
                }
            }
            dolfin::log (dolfin::DBG, 
                         "Removed %d entries from links map", 
                         erasedCount);
            
            dolfin::end ();
        }
    }



    void AbstractEquationSystem::reorderProblems (const std::vector<std::string>& solveOrder)
    {
        dolfin::log (dolfin::DBG, "Setting problems order...");
        solveOrder_ = solveOrder;
    }



    void AbstractEquationSystem::addLink (const std::string& linkFrom, 
                                          const std::string& linkedCoefficientName,
                                          const std::string& linkedCoefficientType, 
                                          const std::string& linkTo,
                                          const bool& forceRelinking)
    {
        dolfin::begin (dolfin::DBG, 
                       "Setting up link (%s, %s, %s) -> (%s, all solution components)...",
                       linkFrom.c_str (),
                       linkedCoefficientName.c_str (),
                       linkedCoefficientType.c_str (),
                       linkTo.c_str ());
        
        // create pair containing the link information passed as input arguments.
        // auto keyword used in place of std::pair <std::tuple <std::string, std::string, std::string>, std::string> 
        // to enhance readability
        auto link = std::make_pair (std::make_tuple (linkFrom, linkedCoefficientName, linkedCoefficientType), 
                                    std::make_pair (linkTo, -1));

        // search for map key in problemsLinks_. 
        // remember that the key (i.e. link.first) is an std::tuple<std::string, std::string, std::string>
        // auto keyword used to enhance readability in place of 
        // std::map <std::tuple <std::string, std::string, std::string>, std::string>::iterator 
        auto linkPosition = problemsLinks_.find (link.first);

        if (linkPosition == problemsLinks_.end ()) // if key not found in map, insert link
        {
            dolfin::log (dolfin::DBG, "Inserting link in links map...");
            problemsLinks_.insert (link);
            
            // perform linking
            linkProblems (link);
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
                                "all solution components)" : 
                                "component " + std::to_string (std::get<1> (linkPosition->second)) + ")")
                << dolfin::endl;

            problemsLinks_.erase (linkPosition);

            dolfin::cout << "and inserting link: " << dolfin::endl;
            dolfin::cout << "\t(" 
                << std::get<0> (link.first)
                << ", " 
                << std::get<1> (link.first)
                << ", " 
                << std::get<2> (link.first)
                << ") -> (" 
                << std::get<0> (link.second)
                << ", all solution components)" 
                << dolfin::endl;

            problemsLinks_.insert (link);
            
            // perform linking
            linkProblems (link);
        }
        else
        {
            dolfin::warning ("link (%s, %s, %s) -> (%s, all solution components) not added. Key is already present in map",
                             (std::get<0> (link.first)).c_str (),
                             (std::get<1> (link.first)).c_str (),
                             (std::get<2> (link.first)).c_str (),
                             (std::get<0> (link.second)).c_str ());
        }
        dolfin::end ();
    }



    void AbstractEquationSystem::addLink (const std::string& linkFrom, 
                                          const std::string& linkedCoefficientName,
                                          const std::string& linkedCoefficientType, 
                                          const std::string& linkTo,
                                          const int& linkToComponent,
                                          const bool& forceRelinking)
    {
        dolfin::begin (dolfin::DBG, 
                       "Setting up link (%s, %s, %s) -> (%s, component %d)...",
                       linkFrom.c_str (),
                       linkedCoefficientName.c_str (),
                       linkedCoefficientType.c_str (),
                       linkTo.c_str (),
                       linkToComponent);
        
        // create pair containing the link information passed as input arguments.
        // auto keyword used in place of 
        // std::pair <std::tuple <std::string, std::string, std::string>, std::pair <std::string, int>> 
        // to enhance readability
        auto link = std::make_pair (std::make_tuple (linkFrom, linkedCoefficientName, linkedCoefficientType), 
                                    std::make_pair (linkTo, linkToComponent));

        // search for map key in problemsLinks_. 
        // remember that the key (i.e. link.first) is an std::tuple<std::string, std::string, std::string>
        // auto keyword used to enhance readability in place of 
        // std::map <std::tuple <std::string, std::string, std::string>, std::pair <std::string, int>>::iterator 
        auto linkPosition = problemsLinks_.find (link.first);

        if (linkPosition == problemsLinks_.end ()) // if key not found in map, insert link
        {
            dolfin::log (dolfin::DBG, "Inserting link in links map...");
            problemsLinks_.insert (link);
            
            // perform linking
            linkProblems (link);
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
                                "all solution components)" : 
                                "component " + std::to_string (std::get<1> (linkPosition->second)) + ")")
                << dolfin::endl;

            problemsLinks_.erase (linkPosition);

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
                << ")"
                << dolfin::endl;

            problemsLinks_.insert (link);
            
            // perform linking
            linkProblems (link);
        }
        else
        {
            dolfin::warning ("link (%s, %s, %s) -> (%s, component %d) not added. Key is already present in map",
                             (std::get<0> (link.first)).c_str (),
                             (std::get<1> (link.first)).c_str (),
                             (std::get<2> (link.first)).c_str (),
                             (std::get<0> (link.second)).c_str (),
                              std::get<1> (link.second));
        }
        dolfin::end ();
    }
            


    const dcp::AbstractProblem& AbstractEquationSystem::operator[] (const std::string& name) const
    {
        auto problemIterator = storedProblems_.find (name);
        if (problemIterator == storedProblems_.end ())
        {
            dolfin::dolfin_error ("dcp: AbstractEquationSystem.cpp",
                                  "operator[]", 
                                  "Problem \"%s\" not found in stored problems map", 
                                  name.c_str ());
        }
        return *(problemIterator->second);
    }



    dcp::AbstractProblem& AbstractEquationSystem::operator[] (const std::string& name)
    {
        auto problemIterator = storedProblems_.find (name);
        if (problemIterator == storedProblems_.end ())
        {
            dolfin::dolfin_error ("dcp: AbstractEquationSystem.cpp",
                                  "operator[]", 
                                  "Problem \"%s\" not found in stored problems map", 
                                  name.c_str ());
        }
        return *(problemIterator->second);
    }



    const dcp::AbstractProblem& AbstractEquationSystem::operator[] (const std::size_t& position) const
    {
        if (position >= solveOrder_.size ())
        {
            dolfin::dolfin_error ("dcp: AbstractEquationSystem.cpp",
                                  "operator[]",
                                  "Input value \"%d\" is greater than problems vector size",
                                  position);
        }
        return this->operator[] (solveOrder_ [position]);
    }



    dcp::AbstractProblem& AbstractEquationSystem::operator[] (const std::size_t& position)
    {
        if (position >= solveOrder_.size ())
        {
            dolfin::dolfin_error ("dcp: AbstractEquationSystem.cpp",
                                  "operator[]",
                                  "Input value \"%d\" is greater than problems vector size",
                                  position);
        }
        return this->operator[] (solveOrder_ [position]);
    }



    void AbstractEquationSystem::print ()
    {
        dolfin::cout << "Problems solve order:" << dolfin::endl;
        for (auto i : solveOrder_)
        {
            dolfin::cout << "\t" << i << dolfin::endl; 
        }
        dolfin::cout << dolfin::endl;

        dolfin::cout << "Problems links:" << dolfin::endl;
        for (auto &i : problemsLinks_)
        {
            dolfin::cout << "("
                << std::get<0> (i.first)
                << ", " 
                << std::get<1> (i.first)
                << ", " 
                << std::get<2> (i.first)
                << ") -> (" 
                << std::get<0> (i.second)
                << ", "
                << std::string (std::get<1> (i.second) == -1 ? 
                                "all solution components)" : 
                                "component " + std::to_string (std::get<1> (i.second)) + ")")
                << dolfin::endl;
        }
    }



    void AbstractEquationSystem::solve (const bool& forceRelinking)
    {
        // this function iterates over solveOrder_ and calls solve (problemName) for each problem, thus delegating
        // to the latter function the task of performing the actual parameters setting and solving
        dolfin::begin ("Solving problems...");
        
        for (auto problem : solveOrder_)
        {
            solve (problem, forceRelinking);
        }
        
        dolfin::end ();
    }



    void AbstractEquationSystem::solve (const std::string& problemName, const bool& forceRelinking)
    {
        dolfin::begin ("Solving problem \"%s\"...", problemName.c_str ());

        // get problem with given name from map. Variable problemIterator will be a
        // std::map <std::string, std::unique_ptr <dcp::AbstractProblem>::iterator
        dolfin::log (dolfin::DBG, 
                     "Looking for problem \"%s\" in problems map...", 
                     problemName.c_str ());
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
        problem.solve ();
        
        dolfin::end ();
    }



    void AbstractEquationSystem::solve (const char* problemName, const bool& forceRelinking)
    {
        solve (std::string (problemName), forceRelinking);
    }



    const dolfin::Function& AbstractEquationSystem::solution (const std::string& problemName) const
    {
        auto problemIterator = storedProblems_.find (problemName);
        if (problemIterator == storedProblems_.end ())
        {
            dolfin::dolfin_error ("dcp: AbstractEquationSystem.cpp",
                                  "solution",
                                  "Problem \"%s\" not found in stored problems map", 
                                  problemName.c_str ());
        }
        return (problemIterator->second)->solution ();
    }
    

    
    /******************* PROTECTED METHODS *******************/
    void AbstractEquationSystem::linkProblems (const std::pair <
                                                                std::tuple <std::string, std::string, std::string>, 
                                                                std::pair  <std::string, int>
                                                               >& link)
    {
        if (std::get<1> (link.second) == -1)
        {
            dolfin::begin (dolfin::DBG, 
                           "Considering link: (%s, %s, %s) -> (%s, all solution componentes)...",
                           (std::get<0> (link.first)).c_str (),
                           (std::get<1> (link.first)).c_str (),
                           (std::get<2> (link.first)).c_str (),
                           (std::get<0> (link.second)).c_str ());
        }
        else
        {
            dolfin::begin (dolfin::DBG, 
                           "Considering link: (%s, %s, %s) -> (%s, component %d)...",
                           (std::get<0> (link.first)).c_str (),
                           (std::get<1> (link.first)).c_str (),
                           (std::get<2> (link.first)).c_str (),
                           (std::get<0> (link.second)).c_str (),
                            std::get<1> (link.second));
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
            dolfin::warning ("Cannot link problem \"%s\" to problem \"%s\". No such problem found in stored problems map",
                             (std::get<0> (link.first)).c_str (),
                             (std::get<0> (link.second)).c_str ());
            dolfin::end ();
            return;
        }

        dcp::AbstractProblem& targetProblem = *(targetProblemIterator->second);

        if (std::get<1> (link.second) == -1)
        {
            dolfin::log (dolfin::DBG, 
                         "Linking coefficient \"%s\" of type \"%s\" of problem \"%s\" to solution of problem \"%s\"...",
                         (std::get<1> (link.first)).c_str (),
                         (std::get<2> (link.first)).c_str (),
                         (std::get<0> (link.first)).c_str (),
                         (std::get<0> (link.second)).c_str ());

            problem.setCoefficient (std::get<2> (link.first), 
                                    dolfin::reference_to_no_delete_pointer (targetProblem.solution ()),
                                    std::get<1> (link.first));
        }
        else
        {
            dolfin::log (dolfin::DBG, 
                         "Linking coefficient \"%s\" of type \"%s\" of problem \"%s\" to component %d of solution of problem \"%s\"...",
                         (std::get<1> (link.first)).c_str (),
                         (std::get<2> (link.first)).c_str (),
                         (std::get<0> (link.first)).c_str (),
                          std::get<1> (link.second),
                         (std::get<0> (link.second)).c_str ());

            int component = std::get<1> (link.second);
            problem.setCoefficient (std::get<2> (link.first), 
                                    dolfin::reference_to_no_delete_pointer (targetProblem.solution ()[component]),
                                    std::get<1> (link.first));

        }
        
        dolfin::end ();
    }
}
