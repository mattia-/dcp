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
        solveType_ ("default"),
        solutionType_ ("default")
    { 
        dolfin::begin (dolfin::DBG, "Building GenericEquationSystem...");
        
        dolfin::log (dolfin::DBG, "Setting up parameters...");
        parameters.add ("subiterations_tolerance", 1e-6);
        parameters.add ("subiterations_maximum_iterations", 10);
            
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
    void GenericEquationSystem::addProblem (const std::string& problemName, 
                                            dcp::GenericProblem& problem)
    {
        if (problemName.empty ())
        {
            dolfin::dolfin_error ("dcp: GenericEquationSystem.cpp",
                                  "addProblem",
                                  "Empty problem names not allowed");
        }

        dolfin::begin (dolfin::DBG, 
                       "Inserting problem \"%s\" in equation system...", 
                       problemName.c_str ());
         
        // insert problem into storedProblems_ 
        dolfin::log (dolfin::DBG, 
                     "Inserting problem in problems map with name \"%s\"...", 
                     problemName.c_str ());
        auto result = storedProblems_.insert 
            (std::make_pair (problemName, std::shared_ptr<dcp::GenericProblem> (problem.clone ())));
        
        // if problem was not inserted in map, issue a warning; else, add it to solveOrder_
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
    
    

    void GenericEquationSystem::addProblem (const std::string& problemName, 
                                            const std::shared_ptr<dcp::GenericProblem> problem)
    {
        if (problemName.empty ())
        {
            dolfin::dolfin_error ("dcp: GenericEquationSystem.cpp",
                                  "addProblem",
                                  "Empty problem names not allowed");
        }

        dolfin::begin (dolfin::DBG, 
                       "Inserting problem \"%s\" in equation system...", 
                       problemName.c_str ());
        
        // insert problem into storedProblems_ 
        dolfin::log (dolfin::DBG, 
                     "Inserting problem in problems map with name \"%s\"...", 
                     problemName.c_str ());
        auto result = storedProblems_.insert (std::make_pair (problemName, problem));
        
        // if problem was not inserted in map, issue a warning; else, add it to vector solveOrder_
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



    void GenericEquationSystem::removeProblem (const std::string& problemName)
    {
        dolfin::begin (dolfin::DBG, 
                       "Removing problem \"%s\" from equation system...", 
                       problemName.c_str ());
        
        // 1) 
        // delete problem from storedProblems_
        dolfin::log (dolfin::DBG, 
                     "Removing problem \"%s\" from problems map...", 
                     problemName.c_str ());
        auto result = storedProblems_.erase (problemName);
        
        // if problem was not found in map, issue a warning; else, erase it also from solveOrder_
        // remember that erase returns the number of elements removed, which in the case of a map is at most 1
        if (result < 1) 
        {
            dolfin::warning ("Problem \"%s\" was not found in equation system. Maybe you used a wrong name?",
                             problemName.c_str ());
        }
        
        // 2)
        // remove problemName from solveOrder_. Remember that we are sure that such name is 
        // unique in vector, since it was either inserted either by the function addProblem 
        // (and if two problems with the same name are inserted, the second one does not get 
        // added either to the map or to the vector) or by the function reorderProblems, 
        // which checks the input vector for double entries when called).
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
    }



    void GenericEquationSystem::reorderProblems (const std::vector<std::string>& solveOrder)
    {
        dolfin::log (dolfin::DBG, "Setting problems order...");

        // look for empty names
        auto emptyName = std::find (solveOrder.begin (), solveOrder.end (), std::string ());
        if (emptyName != solveOrder.end ())
        {
            dolfin::warning ("Input vector to reorderProblems() contains empty names. No reordering performed");
            return;
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
            return;
        }
        
        solveOrder_ = solveOrder;
    }



    void GenericEquationSystem::addLink (const std::string& linkFrom, 
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
            
            // do not perform linking yet. It will be performed once solve is called
            // linkProblems (link);
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
            
            // do not perform linking yet. It will be performed once solve is called
            // linkProblems (link);
        }
        else
        {
            dolfin::warning ("Link (%s, %s, %s) -> (%s, all solution components) not added. Key is already present in map",
                             linkFrom.c_str (),
                             linkedCoefficientName.c_str (),
                             linkedCoefficientType.c_str (),
                             linkTo.c_str ());
        }
        dolfin::end ();
    }



    void GenericEquationSystem::addLink (const std::string& linkFrom, 
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
            
            // do not perform linking yet. It will be performed once solve is called
            // linkProblems (link);
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
            
            // do not perform linking yet. It will be performed once solve is called
            // linkProblems (link);
        }
        else
        {
            dolfin::warning ("Link (%s, %s, %s) -> (%s, component %d) not added. Key is already present in map",
                             linkFrom.c_str (),
                             linkedCoefficientName.c_str (),
                             linkedCoefficientType.c_str (),
                             linkTo.c_str (),
                             linkToComponent);
        }
        dolfin::end ();
    }
    


    bool GenericEquationSystem::removeLink (const std::string& linkFrom, 
                                            const std::string& linkedCoefficientName, 
                                            const std::string& linkedCoefficientType)
    {
        auto linkKey = std::make_tuple (linkFrom, linkedCoefficientName, linkedCoefficientType);
        return (problemsLinks_.erase (linkKey)) == 1 ? true : false;
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
    void GenericEquationSystem::linkProblems (const dcp::GenericEquationSystem::Link& link)
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

        dcp::GenericProblem& problem = *(problemIterator->second);

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

        dcp::GenericProblem& targetProblem = *(targetProblemIterator->second);

        if (std::get<1> (link.second) == -1)
        {
            dolfin::log (dolfin::DBG, 
                         "Linking coefficient \"%s\" of type \"%s\" of problem \"%s\" to solution of problem \"%s\"...",
                         (std::get<1> (link.first)).c_str (),
                         (std::get<2> (link.first)).c_str (),
                         (std::get<0> (link.first)).c_str (),
                         (std::get<0> (link.second)).c_str ());

            // we use solutionType_. Remember that in subiterate() it was set to "stashed", so when relinking takes 
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
            
            // we use solutionType_. Remember that in subiterate() it was set to "stashed", so when relinking takes 
            // place during subiterations the stashed solution will be used
            problem.setCoefficient (std::get<2> (link.first), 
                                    dolfin::reference_to_no_delete_pointer (targetProblem.solution (solutionType_)[component]),
                                    std::get<1> (link.first));

        }
        
        dolfin::end ();
    }
    


    void GenericEquationSystem::subiterate (std::vector<std::string>::const_iterator subiterationsBegin,
                                            std::vector<std::string>::const_iterator subiterationsEnd)
    {
        dolfin::begin (dolfin::PROGRESS, 
                       "===== Subiterating on problems [%s - %s) =====", 
                       ("\"" + *subiterationsBegin + "\"").c_str (),
                       subiterationsEnd == solveOrder_.end () ?  "end" : 
                                                                 ("\"" + *subiterationsEnd + "\"").c_str ());
        
        // get backups for solveType_ and solutionType_
        std::string solveTypeBackup = solveType_;
        std::string solutionTypeBackup = solutionType_;
        
        // set solveType_ so that during subiterations solve of type "stash" is called
        solveType_ = "stash";
        
        // solve all the problems the first time. Remember that subiterationsBegin and subiterationsEnd are
        // iterators on solveOrder_
        dolfin::begin (dolfin::PROGRESS, "Solving all problems once...");
        for (auto problemName = subiterationsBegin; problemName != subiterationsEnd; problemName++)
        {
            solve (*problemName);
        }
        dolfin::end (); // Solving all problems once
            
        // set solution type so that from now on links are performed on stashed solutions
        solutionType_ = "stashed";
        
        // vector to contain auxiliary solutions needed to compute the norm of the increment. Size is equal to
        // the number of problems on which we are subiterating
        dolfin::log (dolfin::DBG, "Setting up subiterations loop variables...");
        std::vector<dolfin::Function> oldSolutions;
        for (auto problemName = subiterationsBegin; problemName != subiterationsEnd; problemName++)
        {
            oldSolutions.push_back (solution (*problemName, solutionType_));
        }
        
        // loop parameters
        double tolerance = parameters ["subiterations_tolerance"];
        int maxIterations = parameters ["subiterations_maximum_iterations"];
        int iteration = 0;
        double sumOfNorms = tolerance + 1;
        dcp::DotProduct dotProduct;
        
        // loop until convergence, that is until the sum of the norms of the increment divided by the norm of the
        // current solution is below tolerance or until maximum number of iterations is reached
        dolfin::begin (dolfin::PROGRESS, "***** SUBITERATIONS LOOP *****");
        while (sumOfNorms >= tolerance && iteration < maxIterations)
        {
            iteration++;
            dolfin::begin (dolfin::PROGRESS, "===== Subiteration %d =====", iteration);
            
            sumOfNorms = 0;
            int counter = 0;
            for (auto problemName = subiterationsBegin; problemName != subiterationsEnd; problemName++)
            {
                solve (*problemName);
                
                dolfin::Function increment ((this -> operator[] (*problemName)).functionSpace ());
                increment = solution (*problemName, solutionType_) - oldSolutions[counter];
                
                oldSolutions[counter] = solution (*problemName, solutionType_);
                
                // use DOLFIN_EPS in division in case norm of the current solution is 0
                sumOfNorms += dotProduct.norm (increment) / (dotProduct.norm (oldSolutions[counter]) + DOLFIN_EPS);
                
                counter++;
            }
            
            dolfin::log (dolfin::PROGRESS, "Sum of relative increment norms: %f", sumOfNorms);
            
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
}
