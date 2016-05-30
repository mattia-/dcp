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

#include <dcp/problems/GenericEquationSystem.h>
#include <dcp/utils/DotProduct.h>
#include <utility>
#include <iterator>
#include <tuple>
#include <dolfin.h>
#include <algorithm>
#include <dolfin/log/dolfin_log.h>

namespace dcp
{
    /******************* CONSTRUCTORS ******************/
    GenericEquationSystem::GenericEquationSystem () : 
        storedProblems_ (),
        solveOrder_ (),
        problemsLinks_ (),
        subiterationsRange_ (),
        dotProducts_ (),
        initialGuessesSetters_ (),
        initialGuessesSettersLinks_ (),
        solveType_ ("default"),
        solutionType_ ("default")
    { 
        dolfin::begin (dolfin::DBG, "Building GenericEquationSystem...");
        
        dolfin::log (dolfin::DBG, "Setting up parameters...");
        parameters.add ("subiterations_tolerance", 1e-6);
        parameters.add ("subiterations_maximum_iterations", 10);
        parameters.add ("plot_subiteration_solutions", false);
        parameters.add ("write_subiteration_solutions_to_file", false);
        parameters.add (dolfin::Parameters ("subiterations_blacklist"));
        parameters.add (dolfin::Parameters ("subiterations_whitelist"));
            
        dolfin::end ();
        
        dolfin::log (dolfin::DBG, "GenericEquationSystem object created");
    }

    

    /******************** GETTERS *********************/
    const std::size_t GenericEquationSystem::size () const
    {
        return storedProblems_.size ();
    }



    const std::vector<std::string>& GenericEquationSystem::problemsNames () const
    {
        return solveOrder_;
    }
    


    const std::pair<std::string, std::string>& GenericEquationSystem::subiterationsRange () const
    {
        return subiterationsRange_;
    }
        
            

    /******************* METHODS *******************/
    bool GenericEquationSystem::addProblem (const std::string& problemName, 
                                            dcp::GenericProblem& problem)
    {
        return addProblem (problemName, 
                           std::shared_ptr<dcp::GenericProblem> (problem.clone ()));
    }
    
    

    bool GenericEquationSystem::addProblem (const std::string& problemName, 
                                            const std::shared_ptr<dcp::GenericProblem> problem)
    {
        dolfin::begin (dolfin::DBG, 
                       "Inserting problem \"%s\" in system problems map...", 
                       problemName.c_str ());
         
        bool problemWasAdded = addProblemToMap_ (storedProblems_, problemName, problem);

        // if problem was added, add it also to solveOrder_ (if problemType is "system")
        if (problemWasAdded == true)
        {
            dolfin::log (dolfin::DBG, 
                         "Inserting problem in problem-names vector with name \"%s\"...", 
                         problemName.c_str ());
            solveOrder_.emplace_back (problemName);
        }

        dolfin::end ();

        return problemWasAdded;
    }



    bool GenericEquationSystem::addInitialGuessSetter (const std::string& name, 
                                                       dcp::GenericProblem& problem)
    {
        return addInitialGuessSetter (name, 
                                      std::shared_ptr<dcp::GenericProblem> (problem.clone ()));
    }
    
    

    bool GenericEquationSystem::addInitialGuessSetter (const std::string& name, 
                                                       const std::shared_ptr<dcp::GenericProblem> problem)
    {
        dolfin::begin (dolfin::DBG, 
                       "Inserting problem \"%s\" in initial guesses setters map...", 
                       name.c_str ());
         
        bool problemWasAdded = addProblemToMap_ (initialGuessesSetters_, name, problem);

        dolfin::end ();

        return problemWasAdded;
    }



    bool GenericEquationSystem::removeProblem (const std::string& problemName)
    {
        dolfin::begin (dolfin::DBG, 
                       "Removing problem \"%s\" from map ...", 
                       problemName.c_str ());
        
        // 1) 
        // delete problem from storedProblems_
        bool problemWasRemoved = removeProblemFromMap_ (storedProblems_, problemName);
        
        // 2)
        // remove problemName from solveOrder_. Remember that we are sure that such name is 
        // unique in vector, since it was either inserted either by the function addProblem (and if two problems with
        // the same name are inserted, the second one does not get added either to the map or to the vector) or by the
        // function reorderProblems, which checks the input vector for double entries when called).
        dolfin::log (dolfin::DBG, "Removing problem \"%s\" from problem-names vector...", problemName.c_str ());
        int erasedCount = 0;

        // std::find will return an iterator to the element if found and vector::end () if not found
        auto problemPosition = find (solveOrder_.begin (), solveOrder_.end (), problemName);
        if (problemPosition != solveOrder_.end ())
        {
            ++erasedCount;
            solveOrder_.erase (problemPosition);
        }
        dolfin::log (dolfin::DBG, "Removed %d entries from problem-names vector", erasedCount);

        // 3)
        // remove problemName from problemsLinks_.
        // Remember (c++ standard):
        // "When erasing from a map, iterators, pointers and references referring to elements removed by the 
        // function are invalidated. All other iterators, pointers and references keep their validity."
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
                
                // Increment iterator to point to next element. This is performed before erasing element, 
                // so that increment is still valid
                ++linksIterator;

                problemsLinks_.erase (auxIterator);
                ++erasedCount;
            }
            else
            {
                ++linksIterator;
            }
        }
        dolfin::log (dolfin::DBG, "Removed %d entries from links map", erasedCount);

        dolfin::end ();

        return problemWasRemoved;
    }



    bool GenericEquationSystem::removeInitialGuessSetter (const std::string& name)
    {
        dolfin::begin (dolfin::DBG, 
                       "Removing problem \"%s\" from map...", 
                       name.c_str ());
        
        // 1) 
        // delete problem from initialGuessesSetters_
        bool problemWasRemoved = removeProblemFromMap_ (initialGuessesSetters_, name);
        
        // 2)
        // remove problemName from initialGuessesSettersLinks_
        // Remember (c++ standard):
        // "When erasing from a map, iterators, pointers and references referring to elements removed by the 
        // function are invalidated. All other iterators, pointers and references keep their validity."
        dolfin::log (dolfin::DBG, 
                     "Removing every occurrence of problem \"%s\" from links map...", 
                     name.c_str ());

        int erasedCount = 0;
        auto linksIterator = initialGuessesSettersLinks_.begin (); // iterator pointing to the first element of the set
        while (linksIterator != initialGuessesSettersLinks_.end ())
        {
            // delete element if problemName appears either as first or as fourth string in the map
            if (std::get<0> (linksIterator->first) == name 
                || 
                std::get<0> (linksIterator->second) == name)
            {
                auto auxIterator = linksIterator; // this will be used for the call to function erase
                
                // Increment iterator to point to next element. This is performed before erasing element, 
                // so that increment is still valid
                ++linksIterator;

                initialGuessesSettersLinks_.erase (auxIterator);
                ++erasedCount;
            }
            else
            {
                ++linksIterator;
            }
        }
        dolfin::log (dolfin::DBG, "Removed %d entries from links map", erasedCount);

        dolfin::end ();

        return problemWasRemoved;
    }



    bool GenericEquationSystem::reorderProblems (const std::vector<std::string>& solveOrder)
    {
        dolfin::log (dolfin::DBG, "Setting problems order...");

        // look for empty names
        auto emptyName = std::find (solveOrder.begin (), solveOrder.end (), std::string ());
        if (emptyName != solveOrder.end ())
        {
            dolfin::warning ("Input vector to reorderProblems() contains empty names. No reordering performed");
            return false;
        }

        // check for duplicates. Basically, sort the vector and delete double entries, then 
        // compare the sizes of this vector with the size of the input vector
        std::vector<std::string> orderedInputVector (solveOrder.begin(), solveOrder.end ());
        std::sort (orderedInputVector.begin (), orderedInputVector.end ());
        
        std::vector<std::string>::iterator it;
        it = std::unique (orderedInputVector.begin(), orderedInputVector.end());

        orderedInputVector.resize (std::distance (orderedInputVector.begin(), it));

        if (solveOrder.size () != orderedInputVector.size ())
        {
            dolfin::warning ("Input vector to reorderProblems() contains duplicates. No reordering performed");
            return false;
        }
        
        solveOrder_ = solveOrder;

        return true;
    }



    bool GenericEquationSystem::addLink (const std::string& linkFrom, 
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
        auto link = std::make_pair (std::make_tuple (linkFrom, linkedCoefficientName, linkedCoefficientType), 
                                    std::make_pair (linkTo, -1));

        bool  linkWasInserted = addLinkToMap_ (problemsLinks_, link, forceRelinking);

        dolfin::end (); // Setting up link

        return linkWasInserted;
    }



    bool GenericEquationSystem::addLink (const std::string& linkFrom, 
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
        
        auto link = std::make_pair (std::make_tuple (linkFrom, linkedCoefficientName, linkedCoefficientType), 
                                    std::make_pair (linkTo, linkToComponent));

        bool linkWasInserted = addLinkToMap_ (problemsLinks_, link, forceRelinking);

        dolfin::end (); // Setting up link

        return linkWasInserted;
    }
    


    bool GenericEquationSystem::removeLink (const std::string& linkFrom, 
                                            const std::string& linkedCoefficientName, 
                                            const std::string& linkedCoefficientType,
                                            const std::string& linkType)
    {
        auto linkKey = std::make_tuple (linkFrom, linkedCoefficientName, linkedCoefficientType);

        dolfin::log (dolfin::DBG, 
                     "Removing link (%s, %s, %s) from \"%s\" map",
                     linkFrom.c_str (),
                     linkedCoefficientName.c_str (),
                     linkedCoefficientType.c_str (),
                     linkType.c_str ());

        bool linkWasRemoved;
        if (linkType == "system")
        {
            linkWasRemoved = removeLinkFromMap_ (problemsLinks_, linkKey);
        }
        else if (linkType == "initial_guess")
        {
            linkWasRemoved = removeLinkFromMap_ (initialGuessesSettersLinks_, linkKey);
        }
        else
        {
            dolfin::dolfin_error ("dcp: GenericEquationSystem.cpp",
                                  "removeLink",
                                  "Unknown link type: %s", 
                                  linkType.c_str ());
            return false;
        }

        return linkWasRemoved;
    }
            


    bool GenericEquationSystem::addInitialGuessSetterLink (const std::string& linkFrom, 
                                                           const std::string& linkedCoefficientName,
                                                           const std::string& linkedCoefficientType, 
                                                           const std::string& linkTo,
                                                           const bool& forceRelinking)
    {
        dolfin::begin (dolfin::DBG, 
                       "Setting up link (%s, %s, %s) -> (%s, all solution components) (initial guess setter)...",
                       linkFrom.c_str (),
                       linkedCoefficientName.c_str (),
                       linkedCoefficientType.c_str (),
                       linkTo.c_str ());
        
        // create pair containing the link information passed as input arguments.
        auto link = std::make_pair (std::make_tuple (linkFrom, linkedCoefficientName, linkedCoefficientType), 
                                    std::make_pair (linkTo, -1));

        bool  linkWasInserted = addLinkToMap_ (initialGuessesSettersLinks_, link, forceRelinking);

        dolfin::end (); // Setting up link

        return linkWasInserted;
    }



    bool GenericEquationSystem::addInitialGuessSetterLink (const std::string& linkFrom, 
                                                           const std::string& linkedCoefficientName,
                                                           const std::string& linkedCoefficientType, 
                                                           const std::string& linkTo,
                                                           const int& linkToComponent,
                                                           const bool& forceRelinking)
    {
        dolfin::begin (dolfin::DBG, 
                       "Setting up link (%s, %s, %s) -> (%s, component %d) (initial guess setter)...",
                       linkFrom.c_str (),
                       linkedCoefficientName.c_str (),
                       linkedCoefficientType.c_str (),
                       linkTo.c_str (),
                       linkToComponent);
        
        auto link = std::make_pair (std::make_tuple (linkFrom, linkedCoefficientName, linkedCoefficientType), 
                                    std::make_pair (linkTo, linkToComponent));

        bool linkWasInserted = addLinkToMap_ (initialGuessesSettersLinks_, link, forceRelinking);

        dolfin::end (); // Setting up link

        return linkWasInserted;
    }
    


    bool GenericEquationSystem::setDotProduct (const std::string& problemName, const dolfin::Form& dotProductComputer)
    {
        dolfin::begin (dolfin::DBG, "Setting dot product computer for problem \"%s\"", problemName.c_str ());

        dcp::DotProduct dotProduct;
        dotProduct.setDotProductComputer (dotProductComputer);

        auto result = dotProducts_.insert (std::pair<std::string, dcp::DotProduct> (problemName, dotProduct));

        dolfin::end (); // "Adding dot product computer for problem \"%s\"", problemName.c_str ()

        return result.second;
    }



    const dcp::GenericProblem& GenericEquationSystem::operator[] (const std::string& name) const
    {
        auto problemIterator = storedProblems_.find (name);
        if (problemIterator == storedProblems_.end ())
        {
            dolfin::dolfin_error ("dcp: GenericEquationSystem.cpp",
                                  "operator[]", 
                                  "Problem \"%s\" not found in stored problems map", 
                                  name.c_str ());
        }
        return *(problemIterator->second);
    }



    dcp::GenericProblem& GenericEquationSystem::operator[] (const std::string& name)
    {
        auto problemIterator = storedProblems_.find (name);
        if (problemIterator == storedProblems_.end ())
        {
            dolfin::dolfin_error ("dcp: GenericEquationSystem.cpp",
                                  "operator[]", 
                                  "Problem \"%s\" not found in stored problems map", 
                                  name.c_str ());
        }
        return *(problemIterator->second);
    }



    const dcp::GenericProblem& GenericEquationSystem::operator[] (const std::size_t& position) const
    {
        if (position >= solveOrder_.size ())
        {
            dolfin::dolfin_error ("dcp: GenericEquationSystem.cpp",
                                  "operator[]",
                                  "Input value \"%d\" is greater than problems vector size",
                                  position);
        }
        return this->operator[] (solveOrder_ [position]);
    }



    dcp::GenericProblem& GenericEquationSystem::operator[] (const std::size_t& position)
    {
        if (position >= solveOrder_.size ())
        {
            dolfin::dolfin_error ("dcp: GenericEquationSystem.cpp",
                                  "operator[]",
                                  "Input value \"%d\" is greater than problems vector size",
                                  position);
        }
        return this->operator[] (solveOrder_ [position]);
    }



    void GenericEquationSystem::setSubiterationRange (const std::string& first, const std::string& last)
    {
        subiterationsRange_.first = first;
        subiterationsRange_.second = last;
    }
    


    void GenericEquationSystem::print ()
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



    const dolfin::Function& GenericEquationSystem::solution (const std::string& problemName,
                                                             const std::string& solutionType) const
    {
        auto problemIterator = storedProblems_.find (problemName);
        if (problemIterator == storedProblems_.end ())
        {
            dolfin::dolfin_error ("dcp: GenericEquationSystem.cpp",
                                  "solution",
                                  "Problem \"%s\" not found in stored problems map", 
                                  problemName.c_str ());
        }
        return (problemIterator->second)->solution (solutionType);
    }
    

    
    void GenericEquationSystem::plotSolution (const std::string& plotType)
    {
        for (auto problemName = solveOrder_.begin (); problemName != solveOrder_.end (); problemName++)
        {
            (this->operator[] (*problemName)).plotSolution (plotType);
        }
    }
            


    void GenericEquationSystem::writeSolutionToFile (const std::string& writeType)
    {
        for (auto problemName = solveOrder_.begin (); problemName != solveOrder_.end (); problemName++)
        {
            (this->operator[] (*problemName)).writeSolutionToFile (writeType);
        }
    }
            


    /******************* PROTECTED METHODS *******************/
    bool GenericEquationSystem::addProblemToMap_ (std::map<std::string, std::shared_ptr <dcp::GenericProblem>>& map,
                                                  const std::string& problemName, 
                                                  const std::shared_ptr<dcp::GenericProblem> problem)
    {
        if (problemName.empty ())
        {
            dolfin::dolfin_error ("dcp: GenericEquationSystem.cpp",
                                  "addProblemToMap_",
                                  "Empty problem names not allowed");
        }

        dolfin::log (dolfin::DBG, 
                     "Inserting problem in problems map with name \"%s\"...", 
                     problemName.c_str ());

        auto result = map.insert (std::make_pair (problemName, problem));

        // if problem was not inserted in map, issue a warning
        if (result.second == false)
        {
            dolfin::warning ("Problem \"%s\" already exist in equation system", 
                             problemName.c_str ());
        }

        return result.second;
    }



    bool GenericEquationSystem::removeProblemFromMap_ (std::map<std::string, std::shared_ptr <dcp::GenericProblem>>& map,
                                                       const std::string& problemName)
    {
        dolfin::log (dolfin::DBG, 
                     "Removing problem \"%s\" from problems map...", 
                     problemName.c_str ());
        auto result = map.erase (problemName);
        
        // if problem was not found in map, issue a warning; else, erase it also from solveOrder_
        // remember that erase returns the number of elements removed, which in the case of a map is at most 1
        if (result < 1) 
        {
            dolfin::warning ("Problem \"%s\" was not found in equation system. Maybe you used a wrong name?",
                             problemName.c_str ());
            return false;
        }

        return true;
    }



    bool GenericEquationSystem::addLinkToMap_ (std::map<dcp::GenericEquationSystem::LinkKey, 
                                                        dcp::GenericEquationSystem::LinkValue>& map, 
                                               const dcp::GenericEquationSystem::Link& link,
                                               const bool& forceRelinking)
    {
        auto linkPosition = map.find (link.first);

        bool linkWasInserted;

        if (linkPosition == map.end ()) // if key not found in map, insert link
        {
            dolfin::log (dolfin::DBG, "Inserting link in links map...");
            auto result = map.insert (link);
            linkWasInserted = result.second;
            
            // do not perform linking yet. It will be performed once solve is called
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

            map.erase (linkPosition);

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

            auto result = map.insert (link);
            linkWasInserted = result.second;
            
            // do not perform linking yet. It will be performed once solve is called
        }
        else
        {
            dolfin::warning ("Link not added. Key is already present in map");
            return false;
        }

        return linkWasInserted;
    }



    bool GenericEquationSystem::removeLinkFromMap_ (std::map<dcp::GenericEquationSystem::LinkKey, 
                                                             dcp::GenericEquationSystem::LinkValue>& map, 
                                                    const dcp::GenericEquationSystem::LinkKey& linkKey)
    {
        return (map.erase (linkKey)) == 1 ? true : false;
    }



    void GenericEquationSystem::linkProblems_ (const dcp::GenericEquationSystem::Link& link,
                                               std::map<std::string, std::shared_ptr <dcp::GenericProblem>>& problemsMap)
    {
        if (std::get<1> (link.second) == -1)
        {
            dolfin::begin (dolfin::DBG, 
                           "Considering link: (%s, %s, %s) -> (%s, all solution components)...",
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
        
        // check if problem whose coefficient we want to set exists
        dolfin::log (dolfin::DBG, 
                     "Looking for problem \"%s\" in problems map...", 
                     (std::get<0> (link.first)).c_str ());
        auto problemIterator = problemsMap.find (std::get<0> (link.first));

        if (problemIterator == problemsMap.end ())
        {
            dolfin::warning ("Problem \"%s\" not found in problems map", 
                             (std::get<0> (link.first)).c_str ());
            dolfin::end ();
            return;
        }

        dcp::GenericProblem& problem = *(problemIterator->second);

        // check if target problem of the link exists
        dolfin::log (dolfin::DBG, "Looking for link target in problems map...");
        auto targetProblemIterator = storedProblems_.find (std::get<0> (link.second));
        if (targetProblemIterator == storedProblems_.end ())
        {
            dolfin::warning ("Cannot link problem \"%s\" to problem \"%s\". No such problem found in problems map",
                             (std::get<0> (link.first)).c_str (),
                             (std::get<0> (link.second)).c_str ());
            dolfin::end ();
            return;
        }

        dcp::GenericProblem& targetProblem = *(targetProblemIterator->second);

        if (std::get<1> (link.second) == -1)
        {
            dolfin::log (dolfin::DBG, 
                         "Linking coefficient \"%s\" of type \"%s\" of problem \"%s\" to solution of problem \"%s\"...",
                         (std::get<1> (link.first)).c_str (),
                         (std::get<2> (link.first)).c_str (),
                         (std::get<0> (link.first)).c_str (),
                         (std::get<0> (link.second)).c_str ());

            // we use solutionType_. Remember that in subiterate_() it was set to "stashed", so when relinking takes 
            // place during subiterations the stashed solution will be used
            problem.setCoefficient (std::get<2> (link.first), 
                                    dolfin::reference_to_no_delete_pointer (targetProblem.solution (solutionType_)),
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
            
            // we use solutionType_. Remember that in subiterate_() it was set to "stashed", so when relinking takes 
            // place during subiterations the stashed solution will be used
            problem.setCoefficient (std::get<2> (link.first), 
                                    dolfin::reference_to_no_delete_pointer (targetProblem.solution (solutionType_)[component]),
                                    std::get<1> (link.first));

        }
        
        dolfin::end ();
    }
    


    void GenericEquationSystem::subiterate_ (std::vector<std::string>::const_iterator subiterationsBegin,
                                             std::vector<std::string>::const_iterator subiterationsEnd)
    {
        dolfin::begin (dolfin::PROGRESS, 
                       "===== Subiterating on problems [%s - %s) =====", 
                       ("\"" + *subiterationsBegin + "\"").c_str (),
                       subiterationsEnd == solveOrder_.end () ?  "end" : 
                                                                 ("\"" + *subiterationsEnd + "\"").c_str ());
        
        // loop parameters
        double tolerance = parameters ["subiterations_tolerance"];
        int maxIterations = parameters ["subiterations_maximum_iterations"];
        bool plotSubiterationSolutions = parameters["plot_subiteration_solutions"];
        bool writeSubiterationSolutions = parameters["write_subiteration_solutions_to_file"];
        dolfin::log (dolfin::DBG, "Tolerance: %f ", tolerance);
        dolfin::log (dolfin::DBG, "Maximum subiterations number: %d ", maxIterations);
        dolfin::log (dolfin::DBG, "Plot subiterations solutions: %s ", plotSubiterationSolutions ? "true" : "false");
        dolfin::log (dolfin::DBG, 
                     "Write subiterations solutions to file: %s ", writeSubiterationSolutions ? "true" : "false");


        // variables to be used later in the function:
        // norms of the increments of the solutions
        std::vector<double> incrementsNorms;
        // problem names whose solution should be used in the convergence check ordered as the subiteration sequence
        std::vector<std::string> sortedConvergenceCheckProblemNames;
        // auxiliary solutions needed to compute the norm of the increment
        std::vector<dolfin::Function> oldSolutions;


        // set initial guesses
        dolfin::begin (dolfin::PROGRESS, "Setting initial guesses...");
        setSubiterationInitialGuesses_ (subiterationsBegin, 
                                        subiterationsEnd,
                                        plotSubiterationSolutions, 
                                        writeSubiterationSolutions);
        dolfin::end (); // Setting initial guesses
            

        // get backups for solveType_ and solutionType_
        std::string solveTypeBackup = solveType_;
        std::string solutionTypeBackup = solutionType_;

        // set solveType_ and solutionType_ so that during subiterations stashed solutions are used
        solveType_ = "stash";
        solutionType_ = "stashed";
        

        // solve all the problems the first time
        dolfin::begin (dolfin::PROGRESS, "Solving all problems once...");
        solveSubiterationProblems_ (subiterationsBegin, 
                                    subiterationsEnd,
                                    0,
                                    plotSubiterationSolutions, 
                                    writeSubiterationSolutions,
                                    incrementsNorms,
                                    sortedConvergenceCheckProblemNames,
                                    oldSolutions,
                                    false);
        dolfin::end (); // Solving all problems once
            

        dolfin::begin (dolfin::DBG, "Setting up subiterations loop variables...");
        setSubiterationLoopVariables_ (subiterationsBegin,
                                       subiterationsEnd, 
                                       sortedConvergenceCheckProblemNames, 
                                       oldSolutions);

        int iteration = 0;
        double maxIncrementNorm = tolerance + 1;

        dolfin::end (); // Setting up subiterations loop variables
        

        // loop until convergence, that is until the max of the norms of the increment divided by the norm of the
        // current solution is below tolerance or until maximum number of iterations is reached
        dolfin::begin (dolfin::PROGRESS, "***** SUBITERATIONS LOOP *****");
        while (maxIncrementNorm >= tolerance && iteration < maxIterations)
        {
            iteration++;
            dolfin::begin (dolfin::PROGRESS, "===== Subiteration %d =====", iteration);
            
            incrementsNorms.clear ();

            solveSubiterationProblems_ (subiterationsBegin, 
                                        subiterationsEnd,
                                        iteration,
                                        plotSubiterationSolutions, 
                                        writeSubiterationSolutions,
                                        incrementsNorms,
                                        sortedConvergenceCheckProblemNames,
                                        oldSolutions,
                                        true);

            maxIncrementNorm = *(std::max_element (incrementsNorms.begin (), incrementsNorms.end ()));
            
            dolfin::log (dolfin::PROGRESS, "Max norm of relative increment: %f", maxIncrementNorm);

            dolfin::begin (dolfin::DBG, "Relative increments norms are:");
            for (auto i = 0; i < sortedConvergenceCheckProblemNames.size (); ++i)
            {
                dolfin::log (dolfin::DBG, 
                             "Problem \"%s\": norm of relative increment = %f", 
                             sortedConvergenceCheckProblemNames[i].c_str (), 
                             incrementsNorms[i]);
            }
            dolfin::end (); // "Relative increments norms are:"

            dolfin::end (); // ===== Iteration =====
        }
        dolfin::end (); // "SUBITERATIONS LOOP"
        
        if (iteration == maxIterations)
        {
            dolfin::warning ("Maximum number of iterations reached in subiterations loop");
        }
        
        // apply stashed solutions
        dolfin::log (dolfin::DBG, "Applying subiterations solutions to problems...");
        for (auto problemName = subiterationsBegin; problemName != subiterationsEnd; problemName++)
        {
            (this -> operator[] (*problemName)).applyStashedSolution ();
        }
        
        // restore original values for solveType_ and solutionType_
        solveType_ = solveTypeBackup;
        solutionType_ = solutionTypeBackup;
        
        dolfin::end (); // "Subiterating on problems begin - end"
    }



    void GenericEquationSystem::setSubiterationInitialGuesses_ 
        (const std::vector<std::string>::const_iterator subiterationsBegin,
         const std::vector<std::string>::const_iterator subiterationsEnd,
         const bool plotSubiterationSolutions,
         const bool writeSubiterationSolutions)
    {
        // remember that subiterationsBegin and subiterationsEnd are iterators on solveOrder_
        for (auto problemName = subiterationsBegin; problemName != subiterationsEnd; problemName++)
        {
            dolfin::begin (dolfin::PROGRESS, "Problem: \"%s\"", (*problemName).c_str ());

            // look for given problem name in initialGuessesSetters_
            auto problemIterator = initialGuessesSetters_.find (*problemName);
            if (problemIterator == initialGuessesSetters_.end ()) // if not found, use current solution as initial guess
            {
                dolfin::log (dolfin::DBG, "Using current solution as initial guess");
                (*this)[*problemName].stashedSolution_ = (*this)[*problemName].solution_.back ().second;
            }
            else // if found, solve the initial guess setter problem
            {
                dolfin::begin (dolfin::DBG, "Solving initial guess setter problem...");
                
                // get problem with given name from map
                dcp::GenericProblem& problem = *(problemIterator->second);

                // 1)
                // loop over initialGuessesSettersLinks_ to reset all links.
                // Remember it is a map: elements in it are order according to the default lexicographic ordering
                dolfin::begin (dolfin::DBG, "Scanning problems links...");

                auto linksIterator = initialGuessesSettersLinks_.begin ();
                while (linksIterator != initialGuessesSettersLinks_.end () 
                       && 
                       std::get<0>(linksIterator->first) <= *problemName)
                {
                    if (std::get<0> (linksIterator->first) == *problemName)
                    {
                        linkProblems_ (*linksIterator, initialGuessesSetters_);
                    }
                    ++linksIterator;
                }

                dolfin::end (); // Scanning problems links

                // 2)
                // solve problem
                problem.solve ();

                // 3) set computed solution as initial guess
                (*this)[*problemName].stashedSolution_ = problem.solution_.back ().second;

                dolfin::end (); // Solving initial guess setter problem
            }

            // plot and write to file
            if (plotSubiterationSolutions == true)
            {
                dolfin::begin (dolfin::DBG, "Plotting subiteration solution...");

                // set plot name to include subiteration number
                std::string oldPlotTitle = (this -> operator[] (*problemName)).parameters["plot_title"];
                std::string newPlotTitle = oldPlotTitle + " (subiteration initial solution)";
                (this -> operator[] (*problemName)).parameters["plot_title"] = newPlotTitle;

                // actual plotting
                (this -> operator[] (*problemName)).plotSolution ("stashed");

                // reset plot title
                (this -> operator[] (*problemName)).parameters["plot_title"] = oldPlotTitle;
                dolfin::end ();
            }

            if (writeSubiterationSolutions == true)
            {
                dolfin::begin ("Writing subiteration solution to file...");
                (this -> operator[] (*problemName)).writeSolutionToFile ("stashed");
                dolfin::end ();
            }

            dolfin::end (); // Problem %s
        }
    }



    void GenericEquationSystem::solveSubiterationProblems_ 
        (const std::vector<std::string>::const_iterator subiterationsBegin,
         const std::vector<std::string>::const_iterator subiterationsEnd,
         const int& iteration,
         const bool plotSubiterationSolutions,
         const bool writeSubiterationSolutions,
         std::vector<double>& incrementsNorms,
         const std::vector<std::string>& sortedConvergenceCheckProblemNames,
         std::vector<dolfin::Function>& oldSolutions,
         const bool inLoop)
    {
        // counter to loop through the convergence check problems. Only used if inLoop is true
        int problemCounter = 0;

        // remember that subiterationsBegin and subiterationsEnd are iterators on solveOrder_
        for (auto problemName = subiterationsBegin; problemName != subiterationsEnd; problemName++)
        {
            dolfin::begin (dolfin::PROGRESS, "Problem: \"%s\"", problemName->c_str ());
            solve_ (*problemName);
            dolfin::end ();
            
            if (inLoop == true)
            {
                if (problemCounter < sortedConvergenceCheckProblemNames.size()
                    && *problemName == sortedConvergenceCheckProblemNames[problemCounter])
                {
                    dolfin::Function increment ((this -> operator[] (*problemName)).functionSpace ());
                    increment = solution (*problemName, solutionType_) - oldSolutions[problemCounter];

                    oldSolutions[problemCounter] = solution (*problemName, solutionType_);

                    // search for current problem in dotProductComputers_
                    auto position = dotProducts_.find (*problemName);

                    // if found, use the dotProduct in the map. Else, create it and use the default
                    if (position != dotProducts_.end ())
                    {
                        dcp::DotProduct& dotProduct = position->second;

                        // use DOLFIN_EPS in division in case norm of the current solution is 0
                        incrementsNorms.push_back 
                            (dotProduct.norm (increment) / (dotProduct.norm (oldSolutions[problemCounter]) + DOLFIN_EPS));
                    }
                    else
                    {
                        dcp::DotProduct dotProduct;

                        // use DOLFIN_EPS in division in case norm of the current solution is 0
                        incrementsNorms.push_back 
                            (dotProduct.norm (increment) / (dotProduct.norm (oldSolutions[problemCounter]) + DOLFIN_EPS));
                    }

                    problemCounter++;
                }
            }

            if (plotSubiterationSolutions == true)
            {
                dolfin::begin (dolfin::DBG, "Plotting subiteration solution...");

                // set plot name to include subiteration number
                std::string oldPlotTitle = (this -> operator[] (*problemName)).parameters["plot_title"];
                std::string newPlotTitle = oldPlotTitle + " (subiteration iteration " + std::to_string (iteration) + ")";
                (this -> operator[] (*problemName)).parameters["plot_title"] = newPlotTitle;

                // actual plotting
                (this -> operator[] (*problemName)).plotSolution ("stashed");

                // reset plot title
                (this -> operator[] (*problemName)).parameters["plot_title"] = oldPlotTitle;
                dolfin::end ();
            }

            if (writeSubiterationSolutions == true)
            {
                dolfin::begin ("Writing subiteration solution to file...");
                (this -> operator[] (*problemName)).writeSolutionToFile ("stashed");
                dolfin::end ();
            }
        }
    }



    void GenericEquationSystem::setSubiterationLoopVariables_ 
        (const std::vector<std::string>::const_iterator subiterationsBegin,
         const std::vector<std::string>::const_iterator subiterationsEnd,
         std::vector<std::string>& sortedConvergenceCheckProblemNames,
         std::vector<dolfin::Function>& oldSolutions)
        {
            // get whitelisted problem names
            std::vector<std::string> convergenceCheckProblemNames;
            parameters ("subiterations_whitelist").get_parameter_keys (convergenceCheckProblemNames);

            // if convergenceCheckProblemNames is empty, there were no whitelisted problems, so use all the problems in the
            // subiteration range minus the blacklisted ones
            if (convergenceCheckProblemNames.empty ())
            {
                for (auto problemName = subiterationsBegin; problemName != subiterationsEnd; problemName++)
                {
                    if (parameters ("subiterations_blacklist").has_parameter (*problemName) == false)
                    {
                        convergenceCheckProblemNames.push_back (*problemName);
                    }
                }
            }

            // sort problems in convergenceCheckProblemNames using the order in the subiteration sequence
            for (auto problemName = subiterationsBegin; problemName != subiterationsEnd; problemName++)
            {
                auto found = std::find (convergenceCheckProblemNames.begin (),
                                        convergenceCheckProblemNames.end (), 
                                        *problemName)
                             != convergenceCheckProblemNames.end ();

                if (found == true)
                {
                    sortedConvergenceCheckProblemNames.push_back (*problemName);
                }
            }

            dolfin::begin (dolfin::DBG, "Problems considered for convergence check are:");
            for (auto& problemName : sortedConvergenceCheckProblemNames)
            {
                dolfin::log (dolfin::DBG, "- " + problemName);
            }
            dolfin::end (); // Problems considered for convergence check are
        
            // get solutions at the initial timestep to be used for the convergence check
            for (auto i = 0; i < sortedConvergenceCheckProblemNames.size(); ++i)
            {
                oldSolutions.push_back (solution (sortedConvergenceCheckProblemNames[i], solutionType_));
            }
        }
}
