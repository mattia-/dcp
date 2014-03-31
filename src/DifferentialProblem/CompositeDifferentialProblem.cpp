#include <DifferentialProblem/CompositeDifferentialProblem.hpp>
#include <utility>
#include <tuple>
#include <dolfin.h>
#include <algorithm>

namespace control_problem
{
    /******************* CONSTRUCTORS ******************/
    CompositeDifferentialProblem::CompositeDifferentialProblem () : 
        storedProblems_ (),
        solveOrder_ (),
        problemsLinks_ (),
        needsLinksScanning_ (true)
    {  }

    

    /******************* METHODS *******************/
    std::size_t CompositeDifferentialProblem::size ()
    {
        return storedProblems_.size ();
    }



    void CompositeDifferentialProblem::addProblem (const std::string& problemName, 
                                                   AbstractDifferentialProblem& problem)
    {
        dolfin::begin (dolfin::DBG, "Inserting problem \"%s\" in composite differential problem...", problemName.c_str ());
        if (dolfin::get_log_level () > dolfin::DBG)
        {
            dolfin::end ();
        }
        
        // create new problem object
        dolfin::log (dolfin::DBG, "Creating new problem object...");
        std::unique_ptr<control_problem::AbstractDifferentialProblem> clonedProblem (problem.clone ());
        
        // insert problem into storedProblems_ taking ownership
        dolfin::log (dolfin::DBG, "Inserting problem in problems map with name \"%s\"...", problemName.c_str ());
        auto result = storedProblems_.insert (std::make_pair (problemName, std::move (clonedProblem)));
        
        // if problem was not inserted in list, issue a warning; else, add it also to the vector solveOrder_
        if (result.second == false)
        {
            dolfin::warning ("Problem \"%s\" already exist in composite differential problem", problemName.c_str ());
        }
        else
        {
            dolfin::log (dolfin::DBG, "Inserting problem in problem-names vector with name \"%s\"...", problemName.c_str ());
            solveOrder_.emplace_back (problemName);
        }
        if (dolfin::get_log_level () <= dolfin::DBG)
        {
            dolfin::end ();
        }
        
    }
    
    void CompositeDifferentialProblem::addProblem (const std::string& problemName, 
                                                   std::unique_ptr<AbstractDifferentialProblem>& problem)
    {
        dolfin::begin (dolfin::DBG, "Inserting problem \"%s\" in composite differential problem...", problemName.c_str ());
        if (dolfin::get_log_level () > dolfin::DBG)
        {
            dolfin::end ();
        }
        
        // insert problem into storedProblems_ taking ownership
        dolfin::log (dolfin::DBG, "Inserting problem in problems map with name \"%s\"...", problemName.c_str ());
        auto result = storedProblems_.insert (std::make_pair (problemName, std::move (problem)));
        
        // if problem was not inserted in list, issue a warning; else, add it also to the vector solveOrder_
        if (result.second == false)
        {
            dolfin::warning ("Problem \"%s\" already exist in composite differential problem", problemName.c_str ());
        }
        else
        {
            dolfin::log (dolfin::DBG, "Inserting problem in problem-names vector with name \"%s\"...", problemName.c_str ());
            solveOrder_.emplace_back (problemName);
        }
        if (dolfin::get_log_level () <= dolfin::DBG)
        {
            dolfin::end ();
        }
    }



    void CompositeDifferentialProblem::removeProblem (const std::string& problemName)
    {
        dolfin::begin (dolfin::DBG, "Removing problem \"%s\" from composite differential problem...", problemName.c_str ());
        
        if (dolfin::get_log_level () > dolfin::DBG)
        {
            dolfin::end ();
        }
        
        // delete problem from storedProblems_
        dolfin::log (dolfin::DBG, "Removing problem \"%s\" from problems map...", problemName.c_str ());
        auto result = storedProblems_.erase (problemName);
        
        // if problem was not inserted in list, issue a warning; else, add it also to the vector solveOrder_
        // remember that erase returns the number of elements removed, which in the case of a map is at most 1
        if (result < 1) 
        {
            if (dolfin::get_log_level () <= dolfin::DBG)
            {
                dolfin::end ();
            }
            dolfin::warning ("Problem \"%s\" was not removed from composite differential problem. Maybe you used a wrong name?", 
                             problemName.c_str ());
        }
        else
        {
            // 1)
            // remove problemName from solveOrder_. Remember that we cannot be sure that such name is 
            // unique in vector and that problems are not in any kind of order.
            // std::find will return an iterator to the element if found and vector::end () if not found
            dolfin::log (dolfin::DBG, "Removing every occurrence of problem \"%s\" from problem-names vector...", 
                         problemName.c_str ());
            int erasedCount = 0;
            auto problemPosition = find (solveOrder_.begin (), solveOrder_.end (), problemName);
            while (problemPosition != solveOrder_.end ())
            {
                ++erasedCount;
                solveOrder_.erase (problemPosition);
                problemPosition = find (solveOrder_.begin (), solveOrder_.end (), problemName);
            }
            dolfin::log (dolfin::DBG, "Removed %d entries from problem-names vector", erasedCount);
            
            // 2)
            // remove problemName from problemsLinks_.
            // We use an important statement from the c++ standard: 
            // When erasing from a map, iterators, pointers and references referring to elements removed by the 
            // function are invalidated. All other iterators, pointers and references keep their validity.
            dolfin::log (dolfin::DBG, "Removing every occurrence of problem \"%s\" from links map...", 
                         problemName.c_str ());
            erasedCount = 0;
            auto linksIterator = problemsLinks_.begin (); // iterator pointing to the first element of the set
            while (linksIterator != problemsLinks_.end ())
            {
                // delete element if problemName appears either as first or as third element in the tuple
                if (std::get<0> (linksIterator->first) == problemName || linksIterator->second == problemName)
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
            dolfin::log (dolfin::DBG, "Removed %d entries from links map", erasedCount);
            
            needsLinksScanning_ = true;
            
            if (dolfin::get_log_level () <= dolfin::DBG)
            {
                dolfin::end ();
            }
        }
    }



    void CompositeDifferentialProblem::reorderProblems (const std::vector<std::string>& solveOrder)
    {
        dolfin::log (dolfin::DBG, "Setting problems order...");
        solveOrder_ = solveOrder;
    }



    void CompositeDifferentialProblem::linkProblems (const std::string& linkFrom, 
                                                     const std::string& linkedCoefficientName,
                                                     const std::string& linkedCoefficientType, 
                                                     const std::string& linkTo,
                                                     const bool& forceRelinking)
    {
        dolfin::begin (dolfin::DBG, "Setting up link (%s, %s, %s) -> %s...",
                       linkFrom.c_str (),
                       linkedCoefficientName.c_str (),
                       linkedCoefficientType.c_str (),
                       linkTo.c_str ());
        
        if (dolfin::get_log_level () > dolfin::DBG)
        {
            dolfin::end ();
        }
        
        // create pair containing the link information passed as input arguments.
        // auto keyword used in place of std::pair <std::tuple <std::string, std::string, std::string>, std::string> 
        // to enhance readability
        auto link = std::make_pair (std::make_tuple (linkFrom, linkedCoefficientName, linkedCoefficientType), linkTo);

        // search for map key in problemsLinks_. 
        // remember that the key (i.e. link.first) is an std::tuple<std::string, std::string, std::string>
        // auto keyword used to enhance readability in place of 
        // std::map <std::tuple <std::string, std::string, std::string>, std::string>::iterator 
        auto linkPosition = problemsLinks_.find (link.first);

        if (linkPosition == problemsLinks_.end ()) // if key not found in map, insert link
        {
            dolfin::log (dolfin::DBG, "Inserting link in links map...");
            problemsLinks_.insert (link);
            if (dolfin::get_log_level () <= dolfin::DBG)
            {
                dolfin::end ();
            }
            
            needsLinksScanning_ = true;
        }
        else if (forceRelinking == true) // if key found in map but forceRelinking set to true, erase 
        // current link and insert the new one
        {
            dolfin::cout << "In composite control problem: erasing link between problems:" << dolfin::endl;
            dolfin::cout << "\t(" 
                << std::get<0> (linkPosition->first) 
                << ", " 
                << std::get<1> (linkPosition->first) 
                << ", " 
                << std::get<2> (linkPosition->first) 
                << ") -> " 
                << linkPosition->second 
                << dolfin::endl;

            problemsLinks_.erase (linkPosition);

            dolfin::cout << "and inserting link: " << dolfin::endl;
            dolfin::cout << "\t(" 
                << std::get<0> (link.first)
                << ", " 
                << std::get<1> (link.first)
                << ", " 
                << std::get<2> (link.first)
                << ") -> " 
                << link.second 
                << dolfin::endl;

            problemsLinks_.insert (link);
            
            needsLinksScanning_ = true;
            
            if (dolfin::get_log_level () <= dolfin::DBG)
            {
                dolfin::end ();
            }
        }
        else
        {
            dolfin::warning ("In composite control problem: Link (%s, %s, %s) -> %s) not added. Key is already present in map",
                             (std::get<0> (link.first)).c_str (),
                             (std::get<1> (link.first)).c_str (),
                             (std::get<2> (link.first)).c_str (),
                             (link.second).c_str ());
            if (dolfin::get_log_level () <= dolfin::DBG)
            {
                dolfin::end ();
            }
        }
    }



    const control_problem::AbstractDifferentialProblem& 
    CompositeDifferentialProblem::operator[] (const std::string& name) const
    {
        auto problemIterator = storedProblems_.find (name);
        if (problemIterator == storedProblems_.end ())
        {
            dolfin::error ("Problem \"%s\" not found in stored problems map", name.c_str ());
        }
        return *(problemIterator->second);
    }



    control_problem::AbstractDifferentialProblem& 
    CompositeDifferentialProblem::operator[] (const std::string& name)
    {
        auto problemIterator = storedProblems_.find (name);
        if (problemIterator == storedProblems_.end ())
        {
            dolfin::error ("Problem \"%s\" not found in stored problems map", name.c_str ());
        }
        return *(problemIterator->second);
    }



    const control_problem::AbstractDifferentialProblem& 
    CompositeDifferentialProblem::operator[] (const std::size_t& position) const
    {
        if (position >= solveOrder_.size ())
        {
            dolfin::error ("Input value \"%d\" is greater than problems vector size", position);
        }
        return this->operator[] (solveOrder_ [position]);
    }



    control_problem::AbstractDifferentialProblem& 
    CompositeDifferentialProblem::operator[] (const std::size_t& position)
    {
        if (position >= solveOrder_.size ())
        {
            dolfin::error ("Input value \"%d\" is greater than problems vector size", position);
        }
        return this->operator[] (solveOrder_ [position]);
    }



    void CompositeDifferentialProblem::print ()
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
                << ") -> " 
                << i.second 
                << dolfin::endl;
        }
    }



    void CompositeDifferentialProblem::solve (const bool& forceRelinking)
    {
        // this function iterates over solveOrder_ and calls solve (problemName) for each problem, thus delegating
        // to the latter function the task of performing the actual parameters setting and solving
        dolfin::begin (dolfin::DBG, "Solving problems...");
        if (dolfin::get_log_level () > dolfin::DBG)
        {
            dolfin::end ();
        }
        
        for (auto problem : solveOrder_)
        {
            solve (problem, forceRelinking);
        }
        
        // at this point, all links were scanned, so needsLinksScanning_ can be set to false. If no operation on
        // links is performed before the next call to solve(), there is no need to scan links again
        needsLinksScanning_ = false;
        
        if (dolfin::get_log_level () <= dolfin::DBG)
        {
            dolfin::end ();
        }
    }



    void CompositeDifferentialProblem::solve (const std::string& problemName, const bool& forceRelinking)
    {
        dolfin::begin (dolfin::PROGRESS, "Solving problem \"%s\"...", problemName.c_str ());
        if (dolfin::get_log_level () > dolfin::PROGRESS)
        {
            dolfin::end ();
        }
        
        if (forceRelinking == true)
        {
            needsLinksScanning_ = true;
            dolfin::log (dolfin::DBG, "Input variable overrides protected member needsLinksScanning_. Value now: %s",
                         needsLinksScanning_ ? "true" : "false");
        }
        
        // get problem with given name from map. Variable problemIterator will be a
        // std::map <std::string, std::unique_ptr <control_problem::AbstractDifferentialProblem>::iterator
        dolfin::log (dolfin::DBG, "Looking for problem \"%s\" in problems map...", problemName.c_str ());
        auto problemIterator = storedProblems_.find (problemName);

        if (problemIterator == storedProblems_.end ())
        {
            dolfin::warning ("Problem \"%s\" not found in stored problems map. Aborting...", problemName.c_str ());
            return;
        }

        control_problem::AbstractDifferentialProblem& problem = *(problemIterator->second);

        // 1)
        // loop over problemsLinks_. Remember it is a map. Elements in it are order according to the default
        // lexicographical ordering
        if (needsLinksScanning_ == true)
        {
            dolfin::begin (dolfin::PROGRESS, "Scanning problems links...");

            if (dolfin::get_log_level () > dolfin::PROGRESS)
            {
                dolfin::end ();
            }

            auto linksIterator = problemsLinks_.begin ();
            while (linksIterator != problemsLinks_.end () && std::get<0> (linksIterator->first) < problemName)
            {
                dolfin::log (dolfin::DBG, 
                             "Considering link: (%s, %s, %s) -> %s...",
                             (std::get<0> (linksIterator -> first)).c_str (),
                             (std::get<1> (linksIterator -> first)).c_str (),
                             (std::get<2> (linksIterator -> first)).c_str (),
                             (linksIterator -> second).c_str ());

                if (std::get<0> (linksIterator->first) == problemName)
                {
                    // check if target problem of the link exists
                    dolfin::log (dolfin::DBG, "Looking for link target in problems map...");
                    auto targetProblemIterator = storedProblems_.find (linksIterator->second);
                    if (targetProblemIterator == storedProblems_.end ())
                    {
                        dolfin::warning ("Cannot link problem \"%s\" to problem \"%s\". No such problem found in \
                                         stored problems map. Aborting...",
                                         problemName.c_str (),
                                         linksIterator->second.c_str ());
                        return;
                    }

                    control_problem::AbstractDifferentialProblem& targetProblem = *(targetProblemIterator->second);

                    dolfin::log (dolfin::PROGRESS, 
                                 "Linking coefficient \"%s\" of type \"%s\" to solution of problem \"%s\"",
                                 (std::get<1> (linksIterator->first)).c_str (),
                                 (std::get<2> (linksIterator->first)).c_str (),
                                 (linksIterator->second).c_str ());

                    problem.setCoefficient (std::get<2> (linksIterator->first), 
                                            dolfin::reference_to_no_delete_pointer (targetProblem.solution ()),
                                            std::get<1> (linksIterator->first));
                }
                else
                {
                    dolfin::log (dolfin::DBG, "Considered link does not match name of problem being solved. Skipping...");
                }

                ++linksIterator;
            }

            if (dolfin::get_log_level () <= dolfin::PROGRESS)
            {
                dolfin::end ();
            }
        }
        else
        {
            dolfin::log (dolfin::DBG, "No links scanning performed. Protected member needsLinksScanning_ value: %s",
                         needsLinksScanning_ ? "true" : "false");
        }

        // 2)
        // solve problem
        dolfin::begin (dolfin::PROGRESS, "Calling solve method for problem \"%s\"...", problemName.c_str ());
        if (dolfin::get_log_level () > dolfin::PROGRESS)
        {
            dolfin::end ();
        }
        
        problem.solve ();
        
        if (dolfin::get_log_level () <= dolfin::PROGRESS)
        {
            dolfin::end ();
        }
            
        if (dolfin::get_log_level () <= dolfin::PROGRESS)
        {
            dolfin::end ();
        }
    }
    


    void CompositeDifferentialProblem::solve (const char* problemName, const bool& forceRelinking)
    {
        solve (std::string (problemName), forceRelinking);
    }
    

   
    const dolfin::Function& CompositeDifferentialProblem::solution (const std::string& problemName) const
    {
        auto problemIterator = storedProblems_.find (problemName);
        if (problemIterator == storedProblems_.end ())
        {
            dolfin::error ("Problem \"%s\" not found in stored problems map", problemName.c_str ());
        }
        return (problemIterator->second)->solution ();
    }
}
