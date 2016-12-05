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

#include <dcp/problems/TimeDependentEquationSystem.h>
#include <dcp/problems/TimeDependentProblem.h>
#include <utility>
#include <tuple>
#include <dolfin.h>
#include <algorithm>
#include <dolfin/log/dolfin_log.h>

namespace dcp
{
    /******************* CONSTRUCTORS ******************/
    TimeDependentEquationSystem::TimeDependentEquationSystem (const std::shared_ptr<dcp::Time>& time,
                                                              const double& startTime,
                                                              const double& dt,
                                                              const double& endTime) :
        GenericEquationSystem (),
        time_ (time),
        startTime_ (startTime),
        dt_ (dt),
        endTime_ (endTime)
    { 
        storedProblemsSolveType_ = "steady";
        dolfin::log (dolfin::DBG, "TimeDependentEquationSystem object created");
    }

    

    /******************* METHODS *******************/
    bool TimeDependentEquationSystem::addProblem (const std::string& problemName, dcp::GenericProblem& problem)
    {
        // try-catch block to check if problem is actually a dcp::TimeDependentProblem&
        try
        {
            dcp::TimeDependentProblem& castProblem = dynamic_cast<dcp::TimeDependentProblem&> (problem);
            
            // if exception was not thrown, check times before reverting back to base class function
            if (time_ != castProblem.time ())
            {
                dolfin::dolfin_error ("dcp: TimeDependentEquationSystem.cpp",
                                      "addProblem",
                                      "Problem's and system's time objects mismatch when adding problem \"%s\" to the map",
                                      problemName.c_str ());
            }
            
            if (dolfin::near (startTime_, castProblem.startTime ()) == 0)
            {
                dolfin::dolfin_error ("dcp: TimeDependentEquationSystem.cpp",
                                      "addProblem",
                                      "Problem's and system's start times mismatch when adding problem \"%s\" to the map",
                                      problemName.c_str ());
            }
            
            if (dolfin::near (endTime_, castProblem.endTime ()) == 0)
            {
                dolfin::dolfin_error ("dcp: TimeDependentEquationSystem.cpp",
                                      "addProblem",
                                      "Problem's and system's end times mismatch when adding problem \"%s\" to the map",
                                      problemName.c_str ());
            }
            
            if (dolfin::near (dt_, castProblem.dt ()) == 0)
            {
                dolfin::dolfin_error ("dcp: TimeDependentEquationSystem.cpp",
                                      "addProblem",
                                      "Problem's and system's time steps mismatch when adding problem \"%s\" to the map",
                                      problemName.c_str ());
            }
            
            // if we got here, everything is fine, so just call the base class function
            return dcp::GenericEquationSystem::addProblem (problemName, problem);
        }
        catch (std::bad_cast& badCast)
        {
            // else, if exception was thrown, it means that problem was not a TimeDependentProblem&, so issue an error
            dolfin::dolfin_error ("dcp: TimeDependentEquationSystem.cpp",
                                  "addProblem",
                                  "Problem given as input argument is not a dcp::TimeDependentProblem");
            return false; // just to suppress the compilation warning
        }
    }
        


    bool TimeDependentEquationSystem::addProblem (const std::string& problemName, 
                                                  const std::shared_ptr<dcp::GenericProblem> problem)
    {
        std::shared_ptr<dcp::TimeDependentProblem> castProblem = 
            std::dynamic_pointer_cast<dcp::TimeDependentProblem> (problem);
        
        // check if cast was successful
        if (castProblem == nullptr)
        {
            // if castProblem is nullptr, it means that problem was not a TimeDependentProblem&, so issue an error
            dolfin::dolfin_error ("dcp: TimeDependentEquationSystem.cpp",
                                  "addProblem",
                                  "Problem given as input argument is not a dcp::TimeDependentProblem");
            return false; // just to suppress the compilation warning
        }
        else
        {
            // if castProblem was not nullptr, check times before reverting back to base class function
            if (time_ != castProblem -> time ())
            {
                dolfin::dolfin_error ("dcp: TimeDependentEquationSystem.cpp",
                                      "addProblem",
                                      "Problem's and system's time objects mismatch when adding problem \"%s\" to the map",
                                      problemName.c_str ());
            }
            
            if (dolfin::near (startTime_, castProblem->startTime ()) == 0)
            {
                dolfin::dolfin_error ("dcp: TimeDependentEquationSystem.cpp",
                                      "addProblem",
                                      "Problem's and system's start times mismatch when adding problem \"%s\" to the map",
                                      problemName.c_str ());
            }
            
            if (dolfin::near (endTime_, castProblem->endTime ()) == 0)
            {
                dolfin::dolfin_error ("dcp: TimeDependentEquationSystem.cpp",
                                      "addProblem",
                                      "Problem's and system's end times mismatch when adding problem \"%s\" to the map",
                                      problemName.c_str ());
            }
            
            if (dolfin::near (dt_, castProblem->dt ()) == 0)
            {
                dolfin::dolfin_error ("dcp: TimeDependentEquationSystem.cpp",
                                      "addProblem",
                                      "Problem's and system's time steps mismatch when adding problem \"%s\" to the map",
                                      problemName.c_str ());
            }
            
            // if we got here, everything is fine, so just call the base class function
            return dcp::GenericEquationSystem::addProblem (problemName, problem);
        }
    }
    


    void TimeDependentEquationSystem::addLinkToPreviousSolution (const std::string& linkFrom, 
                                                                 const std::string& linkedCoefficientName,
                                                                 const std::string& linkedCoefficientType, 
                                                                 const std::string& linkTo,
                                                                 const int& nStepsBack,
                                                                 const bool& forceRelinking)
    {
        // check if nStepsBack is a valid value
        if (nStepsBack < 1)
        {
            dolfin::dolfin_error ("dcp: TimeDependentEquationSystem.cpp",
                                  "addLinkToPreviousSolution",
                                  "Value of nStepsBack should be at least 1. Use addLink() to link to problems at the current timestep");
        }
        
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
            
            // do not perform linking yet. It will be performed once solve is called
            // linkProblemToPreviousSolution_ (link);
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
            
            // do not perform linking yet. It will be performed once solve is called
            // linkProblemToPreviousSolution_ (link);
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
        // check if nStepsBack is a valid value
        if (nStepsBack < 1)
        {
            dolfin::dolfin_error ("dcp: TimeDependentEquationSystem.cpp",
                                  "addLinkToPreviousSolution",
                                  "Value of nStepsBack should be at least 1. Use addLink() to link to problems at the current timestep");
        }
        
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
            
            // do not perform linking yet. It will be performed once solve is called
            // linkProblemToPreviousSolution_ (link);
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
            
            // do not perform linking yet. It will be performed once solve is called
            // linkProblemToPreviousSolution_ (link);
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
    


    const dcp::TimeDependentProblem& TimeDependentEquationSystem::operator[] (const std::string& problemName) const
    {
        auto problemIterator = storedProblems_.find (problemName);
        if (problemIterator == storedProblems_.end ())
        {
            dolfin::dolfin_error ("dcp: TimeDependentEquationSystem.cpp",
                                  "use operator[]", 
                                  "Problem \"%s\" not found in stored problems map", 
                                  problemName.c_str ());
        }
        return *(std::static_pointer_cast<dcp::TimeDependentProblem> (problemIterator->second));
    }



    dcp::TimeDependentProblem& TimeDependentEquationSystem::operator[] (const std::string& problemName)
    {
        auto problemIterator = storedProblems_.find (problemName);
        if (problemIterator == storedProblems_.end ())
        {
            dolfin::dolfin_error ("dcp: TimeDependentEquationSystem.cpp",
                                  "use operator[]", 
                                  "Problem \"%s\" not found in stored problems map", 
                                  problemName.c_str ());
        }
        return *(std::static_pointer_cast<dcp::TimeDependentProblem> (problemIterator->second));
    }



    const dcp::TimeDependentProblem& TimeDependentEquationSystem::operator[] (const std::size_t& position) const
    {
        if (position >= solveOrder_.size ())
        {
            dolfin::dolfin_error ("dcp: TimeDependentEquationSystem.cpp",
                                  "use operator[]",
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
                                  "use operator[]",
                                  "Input value \"%d\" is greater than problems vector size",
                                  position);
        }
        return this->operator[] (solveOrder_ [position]);
    }

    

    std::shared_ptr<const dcp::Time> TimeDependentEquationSystem::time () const
    {
        return time_;
    }

        

    const double& TimeDependentEquationSystem::startTime () const
    {
        return startTime_;
    }

        

    const double& TimeDependentEquationSystem::dt () const
    {
        return dt_;
    }
            


    const double& TimeDependentEquationSystem::endTime () const
    {
        return endTime_;
    }

        

    void TimeDependentEquationSystem::advanceTime ()
    {
        // set new time value
        for (auto& element : storedProblems_)
        {
            // use static_pointer_cast because if the problem was added to the sytem, it is a dcp::TimeDependentProblem
            // for sure (see add() method)
            const std::shared_ptr<dcp::TimeDependentProblem> problem =
                std::static_pointer_cast<dcp::TimeDependentProblem> (element.second);

            problem->advanceTime ();

        }
    }
    


    void TimeDependentEquationSystem::print ()
    {
        dcp::GenericEquationSystem::print ();

        dolfin::cout << "Problems links to previous solutions:" << dolfin::endl;
        for (auto &i : linksToPreviousSolutions_)
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
                                "all solution components, " : 
                                "component " + std::to_string (std::get<1> (i.second)) + ", ")
                << "time step: current "  << std::get<2> (i.second)
                << dolfin::endl;
        }
    }



    void TimeDependentEquationSystem::solve (const std::string& solveType)
    {
        // this function iterates over solveOrder_ and calls solve (problemName) for each problem, thus delegating
        // to the latter function the task of performing the actual parameters setting and solving.

        // check solveType value
        if (solveType != "default" && solveType != "step")
        {
            dolfin::dolfin_error ("dcp: TimeDependentEquationSystem.cpp", 
                                  "solve",
                                  "Unknown solve type \"%s\" requested",
                                  solveType.c_str ());
        }
        
        dolfin::log (dolfin::DBG, "Selected solve type: %s", solveType.c_str ());

        // find subiterations problems
        auto subiterationsBegin = std::find (solveOrder_.begin (), solveOrder_.end (), subiterationsRange_.first);
        auto subiterationsEnd = std::find (solveOrder_.begin (), solveOrder_.end (), subiterationsRange_.second);

        // call right method depending on solveType
        if (solveType == "default")
        {
            solveLoop_ (subiterationsBegin, subiterationsEnd);
        }
        else if (solveType == "step")
        {
            step_ (subiterationsBegin, subiterationsEnd);
        }
    }



    void TimeDependentEquationSystem::setInitialSolution (const std::string& problemName,
                                                          const dolfin::Function& initialSolution, 
                                                          const unsigned int& stepNumber)
    {
        dolfin::begin (dolfin::DBG, "Setting initial solution...");
        (this->operator[](problemName)).setInitialSolution (initialSolution, stepNumber);
        dolfin::end (); // "Setting initial solution"
    }
    


    void TimeDependentEquationSystem::setInitialSolution (const std::string& problemName,
                                                          const dolfin::Expression& initialSolution, 
                                                          const unsigned int& stepNumber)
    {
        dolfin::begin (dolfin::DBG, "Setting initial solution...");
        (this->operator[](problemName)).setInitialSolution (initialSolution, stepNumber);
        dolfin::end (); // "Setting initial solution"
    }



    void TimeDependentEquationSystem::setInitialSolution (const std::string& problemName)
    {
        std::string storedProblemsSolveTypeBackup = storedProblemsSolveType_;
        storedProblemsSolveType_ = "steady";
        dolfin::begin (dolfin::DBG, "Setting initial solution...");
        solve_ (problemName);
        dolfin::end ();
        storedProblemsSolveType_ = storedProblemsSolveTypeBackup;
    }



    /******************* PROTECTED METHODS *******************/
    void TimeDependentEquationSystem::solve_ (const std::string& problemName)
    {
        // get problem with given name from map
        dcp::TimeDependentProblem& problem = this -> operator[] (problemName);

        // 1)
        // loop over problemsLinks_ to reset all links to take changes to coefficients 
        // into account. Remember it is a map: elements in it are order according to 
        // the default lexicographical ordering
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

        auto previousSolutionsLinksIterator = linksToPreviousSolutions_.begin ();
        while (previousSolutionsLinksIterator != linksToPreviousSolutions_.end () 
               && 
               std::get<0> (previousSolutionsLinksIterator->first) <= problemName)
        {
            if (std::get<0> (previousSolutionsLinksIterator->first) == problemName)
            {
                linkProblemToPreviousSolution_ (*previousSolutionsLinksIterator);
            }
            ++previousSolutionsLinksIterator;
        }

        dolfin::end ();

        // 2)
        // solve problem
        problem.solve (storedProblemsSolveType_);
    }



    void TimeDependentEquationSystem::linkProblemToPreviousSolution_ 
        (const dcp::TimeDependentEquationSystem::PreviousSolutionLink& link)
    {
        if (std::get<1> (link.second) == -1)
        {
        dolfin::begin (dolfin::DBG, 
                           "Considering link: (%s, %s, %s) -> (%s, all solution components, %d time steps back)...",
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

        dcp::GenericProblem& problem = *(problemIterator->second);

        // check if target problem of the link exists
        dolfin::log (dolfin::DBG, "Looking for link target in problems map...");
        auto targetProblemIterator = storedProblems_.find (std::get<0> (link.second));
        if (targetProblemIterator == storedProblems_.end ())
        {
            dolfin::warning ("Cannot link problem \"%s\" to previous solution of problem \"%s\". Target problem not found in problems' map",
                             (std::get<0> (link.first)).c_str (),
                             (std::get<0> (link.second)).c_str ());
            dolfin::end ();
            return;
        }

        // unlike in GenericEquationSystem, we use a reference to TimeDependentProblem (i.e. the derived class).
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
            auto& targetProblemSolutionsVector = targetProblem.solutionsVector ();  
            
            // get target function, by going back from the last element of nStepsBack steps. 
            // NB: we use operator+ to traverse the vector backwards, since rbegin is a REVERSE iterator
            int nStepsBack = std::get<2> (link.second);
            
            // if we are linking a problem (in a subiterations loop) to the solution at a previous problem (in the
            // subiterations loop as well), we must decrease the number of time steps by one before linking because
            // the solution at the current step is not actually stored in solutionsVector. This means that
            // for example the solution at the previous time step is the LAST element in solutionsVector, not the last 
            // but one.
            // So check if they both are between subiterationsRange_ and if that's the case decrease nStepsBack by 1
            auto subiterationsBegin = std::find (solveOrder_.begin (), solveOrder_.end (), subiterationsRange_.first);
            auto subiterationsEnd = std::find (solveOrder_.begin (), solveOrder_.end (), subiterationsRange_.second);
            bool bothInSubiterationsRange = 
                (std::find (subiterationsBegin, subiterationsEnd, problemIterator->first) != solveOrder_.end ())
                && (std::find (subiterationsBegin, subiterationsEnd, targetProblemIterator->first) != solveOrder_.end ());
            if (bothInSubiterationsRange)
            {
                nStepsBack--;
            }
            
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
            auto& targetProblemSolutionsVector = targetProblem.solutionsVector ();  
            
            // get target function, by going back from the last element of nStepsBack steps. 
            // NB: we use operator+ to traverse the vector backwards, since rbegin is a REVERSE iterator
            int nStepsBack = std::get<2> (link.second);
            
            // if we are linking a problem (in a subiterations loop) to the solution at a previous problem (in the
            // subiterations loop as well), we must decrease the number of time steps by one before linking because
            // the solution at the current step is not actually stored in solutionsVector. This means that
            // for example the solution at the previous time step is the LAST element in solutionsVector, not the last 
            // but one.
            // So check if they both are between subiterationsRange_ and if that's the case decrease nStepsBack by 1
            auto subiterationsBegin = std::find (solveOrder_.begin (), solveOrder_.end (), subiterationsRange_.first);
            auto subiterationsEnd = std::find (solveOrder_.begin (), solveOrder_.end (), subiterationsRange_.second);
            bool bothInSubiterationsRange = 
                (std::find (subiterationsBegin, subiterationsEnd, problemIterator->first) != solveOrder_.end ())
                && (std::find (subiterationsBegin, subiterationsEnd, targetProblemIterator->first) != solveOrder_.end ());
            if (bothInSubiterationsRange)
            {
                nStepsBack--;
            }
            
            const dolfin::Function& targetFunction = (targetProblemSolutionsVector.rbegin() + nStepsBack)->second;

            int component = std::get<1> (link.second);
            problem.setCoefficient (std::get<2> (link.first), 
                                    dolfin::reference_to_no_delete_pointer (targetFunction [component]),
                                    std::get<1> (link.first));

        }
        
        dolfin::end ();
    }



    void TimeDependentEquationSystem::solveLoop_ (const std::vector<std::string>::const_iterator subiterationsBegin, 
                                                  const std::vector<std::string>::const_iterator subiterationsEnd)
    {
        // Reserve space for solutions vector of each problem
        for (const auto& element : storedProblems_)
        {
            std::string problemName = element.first;
            dcp::TimeDependentProblem& problem = this -> operator[] (problemName);
            problem.reserve();
        }

        // Solutions loop
        dolfin::begin (dolfin::PROGRESS, "Start solution loop...");
        std::size_t timestep = 0;
        while (isFinished () == 0)
        {
            timestep++;
            dolfin::begin (dolfin::PROGRESS, "===== Timestep %d =====", timestep);
            
            step_ (subiterationsBegin, subiterationsEnd);

            // plot and write to file problems solutions
            for (const auto& problemName : solveOrder_)
            {
                // plot solution
                dcp::TimeDependentProblem& problem = this -> operator[] (problemName);
                int plotInterval = problem.parameters["plot_interval"];

                if (plotInterval > 0 && timestep % plotInterval == 0)
                {
                    problem.plotSolution ("last");
                }

                // write solution to file 
                int writeInterval = problem.parameters ["write_interval"];
                if (writeInterval > 0 && timestep % writeInterval == 0)
                {
                    problem.writeSolutionToFile ("last");
                }
            }

            dolfin::end (); // Timestep %d
        }

        dolfin::end (); // Start solution loop...
    }
        


    void TimeDependentEquationSystem::step_ (const std::vector<std::string>::const_iterator subiterationsBegin, 
                                             const std::vector<std::string>::const_iterator subiterationsEnd)
    {
        advanceTime ();

        dolfin::log (dolfin::PROGRESS, "TIME = %f s", time_ -> value ());

        auto problemName = solveOrder_.cbegin ();
        while (problemName != solveOrder_.cend ())
        {
            if (problemName != subiterationsBegin)
            {
                // solve problem
                dolfin::begin (dolfin::PROGRESS, "Problem: \"%s\"", problemName->c_str ());
                solve_ (*problemName);
                dolfin::end ();

                problemName++;
            }
            else
            {
                // subiterate on problems
                subiterate_ (subiterationsBegin, subiterationsEnd);

                // set problemName to subiterationsEnd, since all problems in subiterations range have already been
                // solved
                problemName = subiterationsEnd;
            }
        }
    }
}
