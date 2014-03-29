#include <DifferentialProblem/CompositeDifferentialProblem.hpp>
#include <utility>
#include <tuple>
#include <dolfin.h>
#include <algorithm>

namespace control_problem
{

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
        
        // if problem was not inserted in list, issue a warning; else, add it also to the vector problemsOrder_
        if (result.second == false)
        {
            dolfin::warning ("Problem \"%s\" already exist in composite differential problem", problemName.c_str ());
        }
        else
        {
            dolfin::log (dolfin::DBG, "Inserting problem in problem-names vector with name \"%s\"...", problemName.c_str ());
            problemsOrder_.emplace_back (problemName);
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
        
        // if problem was not inserted in list, issue a warning; else, add it also to the vector problemsOrder_
        if (result.second == false)
        {
            dolfin::warning ("Problem \"%s\" already exist in composite differential problem", problemName.c_str ());
        }
        else
        {
            dolfin::log (dolfin::DBG, "Inserting problem in problem-names vector with name \"%s\"...", problemName.c_str ());
            problemsOrder_.emplace_back (problemName);
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
        
        // if problem was not inserted in list, issue a warning; else, add it also to the vector problemsOrder_
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
            // remove problemName from problemsOrder_. Remember that we cannot be sure that such name is 
            // unique in vector and that problems are not in any kind of order.
            // std::find will return an iterator to the element if found and vector::end () if not found
            dolfin::log (dolfin::DBG, "Removing every occurrence of problem \"%s\" from problem-names vector...", 
                         problemName.c_str ());
            int erasedCount = 0;
            auto problemPosition = find (problemsOrder_.begin (), problemsOrder_.end (), problemName);
            while (problemPosition != problemsOrder_.end ())
            {
                ++erasedCount;
                problemsOrder_.erase (problemPosition);
                problemPosition = find (problemsOrder_.begin (), problemsOrder_.end (), problemName);
            }
            dolfin::log (dolfin::DBG, "Removed %d entries from problem-names vector", erasedCount);
            
            // 2)
            // remove problemName from linkedProblems_.
            // We use an important statement from the c++ standard: 
            // When erasing from a map, iterators, pointers and references referring to elements removed by the 
            // function are invalidated. All other iterators, pointers and references keep their validity.
            dolfin::log (dolfin::DBG, "Removing every occurrence of problem \"%s\" from links map...", 
                         problemName.c_str ());
            erasedCount = 0;
            auto linksIterator = linkedProblems_.begin (); // iterator pointing to the first element of the set
            while (linksIterator != linkedProblems_.end ())
            {
                // delete element if problemName appears either as first or as third element in the tuple
                if (std::get<0> (linksIterator->first) == problemName || linksIterator->second == problemName)
                {
                    auto auxIterator = linksIterator; // this will be used for the call to function erase
                    ++linksIterator;   // iterator incremented to point to next element. 
                                       // This is performed before erasing element, 
                                       // so that increment is still valid
                    
                    linkedProblems_.erase (auxIterator);
                    ++erasedCount;
                }
                else
                {
                    ++linksIterator;
                }
            }
            dolfin::log (dolfin::DBG, "Removed %d entries from links map", erasedCount);
            
            if (dolfin::get_log_level () <= dolfin::DBG)
            {
                dolfin::end ();
            }
        }
    }



    void CompositeDifferentialProblem::reorderProblems (const std::vector<std::string>& problemsOrder)
    {
        dolfin::log (dolfin::DBG, "Setting problems order...");
        problemsOrder_ = problemsOrder;
    }



    void CompositeDifferentialProblem::linkProblems (const std::string& linkFrom, 
                                                     const std::string& linkedCoefficientName,
                                                     const std::string& linkedCoefficientType, 
                                                     const std::string& linkTo,
                                                     const bool& forceRelinking = false)
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

        // search for map key in linkedProblems_. 
        // remember that the key (i.e. link.first) is an std::tuple<std::string, std::string, std::string>
        // auto keyword used to enhance readability in place of 
        // std::map <std::tuple <std::string, std::string, std::string>, std::string>::iterator 
        auto linkPosition = linkedProblems_.find (link.first);

        if (linkPosition == linkedProblems_.end ()) // if key not found in map, insert link
        {
            dolfin::log (dolfin::DBG, "Inserting link in links map...");
            linkedProblems_.insert (link);
            if (dolfin::get_log_level () <= dolfin::DBG)
            {
                dolfin::end ();
            }
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

            linkedProblems_.erase (linkPosition);

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

            linkedProblems_.insert (link);
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
        if (position >= problemsOrder_.size ())
        {
            dolfin::error ("Input value \"%d\" is greater than problems vector size", position);
        }
        return this->operator[] (problemsOrder_ [position]);
    }



    control_problem::AbstractDifferentialProblem& 
    CompositeDifferentialProblem::operator[] (const std::size_t& position)
    {
        if (position >= problemsOrder_.size ())
        {
            dolfin::error ("Input value \"%d\" is greater than problems vector size", position);
        }
        return this->operator[] (problemsOrder_ [position]);
    }



    void CompositeDifferentialProblem::print ()
    {
        dolfin::cout << "Problems solve order:" << dolfin::endl;
        for (auto i : problemsOrder_)
        {
            dolfin::cout << "\t" << i << dolfin::endl; 
        }
        dolfin::cout << dolfin::endl;

        dolfin::cout << "Problems links:" << dolfin::endl;
        for (auto &i : linkedProblems_)
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



    void CompositeDifferentialProblem::solve ()
    {
        // this function iterates over problemsOrder_ and calls solve (problemName) for each problem, thus delegating
        // to the latter function the task of performing the actual parameters setting and solving
        dolfin::begin (dolfin::DBG, "Solving problems...");
        if (dolfin::get_log_level () > dolfin::DBG)
        {
            dolfin::end ();
        }
        
        for (auto problem : problemsOrder_)
        {
            solve (problem);
        }
        
        if (dolfin::get_log_level () <= dolfin::DBG)
        {
            dolfin::end ();
        }
    }



    void CompositeDifferentialProblem::solve (const std::string& problemName)
    {
        dolfin::begin (dolfin::PROGRESS, "Solving problem \"%s\"...", problemName.c_str ());
        if (dolfin::get_log_level () > dolfin::PROGRESS)
        {
            dolfin::end ();
        }
        
        // get problem with given name from map. Variable problemIterator will be a
        // std::map <std::string, std::unique_ptr <control_problem::AbstractDifferentialProblem>::iterator
        dolfin::log (dolfin::DBG, "Looking for problem \"%s\" in problems map", problemName.c_str ());
        auto problemIterator = storedProblems_.find (problemName);

        if (problemIterator == storedProblems_.end ())
        {
            dolfin::warning ("Problem \"%s\" not found in stored problems map. Aborting...", problemName.c_str ());
            return;
        }

        control_problem::AbstractDifferentialProblem& problem = *(problemIterator->second);

        // 1)
        // loop over linkedProblems_. Remember it is a map. Elements in it are order according to the default
        // lexicographical ordering
        dolfin::begin (dolfin::PROGRESS, "Scanning problems links...");

        if (dolfin::get_log_level () > dolfin::PROGRESS)
        {
            dolfin::end ();
        }
        
        auto linksIterator = linkedProblems_.begin ();
        while (linksIterator != linkedProblems_.end () && std::get<0> (linksIterator->first) < problemName)
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
                                        boost::shared_ptr<dolfin::Function> (new dolfin::Function (targetProblem.solution ())),
                                        std::get<1> (linksIterator->first));
            }
            ++linksIterator;
        }
        
        if (dolfin::get_log_level () <= dolfin::PROGRESS)
        {
            dolfin::end ();
        }

        // 2)
        // solve problem
        dolfin::begin (dolfin::PROGRESS, "Calling solve method for problem \"%s\"", problemName.c_str ());
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
