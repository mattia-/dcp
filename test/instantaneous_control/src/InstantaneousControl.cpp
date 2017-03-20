/* 
 *  Copyright (C) 2017, Ivan Fumagalli, ivan.fumagalli@polimi.it
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

#include "InstantaneousControl.h"// <dcp/differential_problems/InstantaneousControl.h>
#include <dcp/differential_problems/TimeDependentProblem.h>
#include <dcp/differential_problems/MovingTimeDependentProblem.h>
#include <utility>
#include <tuple>
#include <dolfin.h>
#include <algorithm>
#include <dolfin/log/dolfin_log.h>

#include <dcp/differential_problems/LinearProblem.h>
#include "navierstokesAdjoint.h"

extern struct ProblemData problemData;

namespace dcp
{
    /******************* CONSTRUCTORS ******************/
    //InstantaneousControl::InstantaneousControl (std::shared_ptr<dolfin::Form> objectiveFunctional, std::shared_ptr<dolfin::GenericFunction> control) :
    InstantaneousControl::InstantaneousControl (std::shared_ptr<dolfin::Form> objectiveFunctional, std::shared_ptr<dolfin::Constant> control) :
        EquationSystem (),
        objectiveFunctional_ (objectiveFunctional),
        control_ (control)
    { 
        dolfin::log (dolfin::DBG, "InstantaneousControl object created");
    }

    //InstantaneousControl::InstantaneousControl (std::shared_ptr<dolfin::Form> objectiveFunctional, std::shared_ptr<dolfin::Form> projector, std::shared_ptr<dolfin::GenericFunction> control) :
    InstantaneousControl::InstantaneousControl (std::shared_ptr<dolfin::Form> objectiveFunctional, std::shared_ptr<dolfin::Form> projector, std::shared_ptr<dolfin::Constant> control) :
        EquationSystem (),
        objectiveFunctional_ (objectiveFunctional),
        projector_ (projector),
        control_ (control)
    { 
        dolfin::log (dolfin::DBG, "InstantaneousControl object created");
    }

    

    /******************* METHODS *******************/
    void InstantaneousControl::addProblem (const std::string& problemName, dcp::AbstractProblem& problem, bool isTimeDependent)
    {
        this->EquationSystem::addProblem (problemName, problem);

        if (isTimeDependent)
          timeDependentStoredProblems_.insert (std::make_pair (problemName, this->storedProblems_[problemName]));
    }

    void InstantaneousControl::addLinkToPreviousSolution (const std::string& linkFrom, 
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

    bool InstantaneousControl::isFinished ()
    {
        std::size_t nFinished = 0;
        for (auto mapElement : timeDependentStoredProblems_)
        {
            dcp::TimeDependentProblem& tmp = static_cast<dcp::TimeDependentProblem&> (*(mapElement.second));
            nFinished += tmp.isFinished ();
        }
        return (nFinished == timeDependentStoredProblems_.size ());
    }
    
    void InstantaneousControl::addLinkToPreviousSolution (const std::string& linkFrom, 
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
//            linkProblemToPreviousSolution (link);
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

    void InstantaneousControl::solve (const bool& forceRelinking)
    {
        // instrumental variables
        dolfin::Array<double> controlValues (2);
        std::array<double, 2> cP ({0.5*problemData.lx,0});
        dolfin::Array<double> controlPoint (2, cP.data());
        std::string controlName (parameters["control_name"]);
        double alpha (parameters["alpha"]);
        double alphaMin (parameters["alpha_min"]);

        // evaluation of the initial value of the functional
        objectiveFunctional_->set_coefficient ("u", dolfin::reference_to_no_delete_pointer ((*this)["primal"].solution()[0]));
        objectiveFunctional_->set_coefficient (controlName, control_);
        rescale_.setFactor (parameters["penaltyFactor"]);
        //rescale_.setTime (0);
        rescale_.setN (0);
        objectiveFunctional_->set_coefficient ("rescale", dolfin::reference_to_no_delete_pointer(rescale_));
        objFunValues_.push_back (dolfin::assemble (* objectiveFunctional_));

        // file-recording of the initial values
        control_->eval (controlValues, controlPoint);
        std::ofstream file_toBeReset ((problemData.savepath+"objFunBefore_ctrl_alpha_objFun.csv").c_str(), std::ofstream::out);
        file_toBeReset.close();
        file_toBeReset.precision (15);
        file_toBeReset << objFunValues_.back() << ',' << controlValues[1] << ',' << alpha << std::endl;
        file_toBeReset.close();
        this->meshManager().print2csv (
                                (*this)["primal"].solution (), 
                                problemData.savepath+"primal", "0",
                                true );
        this->meshManager().print2csv (
                                (*this)["adjoint"].solution (), 
                                problemData.savepath+"adjoint", "0",
                                true );

        // this function iterates over solveOrder_ and calls solve (problemName) for each problem, thus delegating
        // to the latter function the task of performing the actual parameters setting and solving.
        // The loop is repeated until isFinished() returns true, that is until all problems' time loops are ended
        dolfin::begin ("Solving problems...");
        
        int iterationCounter = 0; 
        while (isFinished () == 0)
        {

            double alpha (parameters["alpha"]);
            double alphaMin (parameters["alpha_min"]);

            iterationCounter++;
            dolfin::log (dolfin::INFO, "=====================================");
            dolfin::log (dolfin::INFO, "TIME DEPENDENT SYSTEM ITERATION: %d", iterationCounter);
            dolfin::begin (dolfin::INFO, "=====================================");

            dolfin::begin ("Solving auxiliary problems...");

            // 1) primal state problem first solution and file-recording
            (*this)["primal"].setCoefficient ("linear_form", control_, controlName);
            solve("primal", forceRelinking, false);
std::cerr << __FILE__ << ' ' << __LINE__ << std::endl;

            this->meshManager().print2csv (
                                    (*this)["primal"].solution (), 
                                    problemData.savepath+"primalBefore",std::to_string(iterationCounter),
                                    true );

            // 2) ALE problem solution (for adjoint problem)
//            solve("ALE", true);//forceRelinking);
            dolfin::Function w (* static_cast<dcp::MovingTimeDependentProblem &> ((*this)["primal"]) .meshManager ().displacement() );
            (*this)["adjoint"].setCoefficient ("bilinear_form", dolfin::reference_to_no_delete_pointer (w), "w");
            const dolfin::Function & primalSol ((*this)["primal"].solution ());
            (*this)["adjoint"].setCoefficient ("linear_form", dolfin::reference_to_no_delete_pointer (primalSol[0]), "u");
            const dolfin::Function & primalOldSol ((static_cast<dcp::TimeDependentProblem &>((*this)["primal"]).solutions ().end() - 1)->second);
            (*this)["adjoint"].setCoefficient ("bilinear_form", dolfin::reference_to_no_delete_pointer (primalOldSol[0]), "u_old");

            // 3) adjoint problem solution and file-recording
            solve("adjoint", false);//forceRelinking);

            this->meshManager().print2csv (
                                    (*this)["adjoint"].solution (), 
                                    problemData.savepath+"adjoint",std::to_string(iterationCounter),
                                    true );

            dolfin::end ();

            // 4) first evaluation of the functional
            objectiveFunctional_->set_coefficient ("u", dolfin::reference_to_no_delete_pointer ((*this)["primal"].solution()[0]));
            objectiveFunctional_->set_coefficient (controlName, control_);
            //rescale_.setTime (static_cast<dcp::TimeDependentProblem&> ((*this)["primal"]).time());
            rescale_.setN (iterationCounter - 1);
            objectiveFunctional_->set_coefficient ("rescale", dolfin::reference_to_no_delete_pointer(rescale_));
            objFunValues_.push_back (dolfin::assemble (* objectiveFunctional_));

            // 5) construction of the control increment and back-tracking
            dolfin::begin("Computing the new control...");
            projector_->set_coefficient ("u", dolfin::reference_to_no_delete_pointer ((*this)["adjoint"].solution ()[0]));
            double projected (dolfin::assemble (* projector_));
            double rescaleValue (rescale_.value());
            control_->eval (controlValues, controlPoint);
            control_.reset (new dolfin::Constant (0, 
                controlValues[1]
                + alpha
                * static_cast<dcp::TimeDependentProblem&> ((*this)["primal"]).dt()
                * ( projected
                    - rescaleValue * problemData.lx * controlValues[1] ) ) );

            alpha = this->chooseIncrementStep (alpha, alphaMin, controlValues, controlPoint, controlName, projected, forceRelinking);

//imposedControl
/*control_.reset (new dolfin::Constant(0, -100e-6 * sin (2*DOLFIN_PI/10 * iterationCounter-1)));
std::ofstream ctrlFile ((problemData.savepath+"ctrlFile.csv").c_str(), std::ofstream::app);
control_->eval (controlValues, controlPoint);
ctrlFile << controlValues[1] << std::endl;*/

            dolfin::end ();

            // 6) solution of the primal problem with the improved control
            dolfin::begin ("Solving the corrected primal problem...");

            (*this)["primal"].setCoefficient ("linear_form", control_, controlName);
            solve("primal", forceRelinking, true);

            this->meshManager().print2csv (
                                    (*this)["primal"].solution (), 
                                    problemData.savepath+"primal",std::to_string(iterationCounter),
                                    true );
            std::ofstream file ((problemData.savepath+"objFunBefore_ctrl_alpha_objFun.csv").c_str(), std::ofstream::app);
            file.precision (15);
            file << objFunValues_.back() << ',' << controlValues[1] << ',' << alpha << ',';
            objectiveFunctional_->set_coefficient ("u", dolfin::reference_to_no_delete_pointer ((*this)["primal"].solution()[0]));
            objectiveFunctional_->set_coefficient (controlName, control_);
            //rescale_.setTime (static_cast<dcp::TimeDependentProblem&> ((*this)["primal"]).time());
            rescale_.setN (iterationCounter - 1);
            objectiveFunctional_->set_coefficient ("rescale", dolfin::reference_to_no_delete_pointer(rescale_));
            objFunValues_.back () = dolfin::assemble (* objectiveFunctional_);
            file << objFunValues_.back() << std::endl;
            file.close();

            dolfin::end ();

            dolfin::log (dolfin::INFO, "");
            dolfin::end ();
        }

        dolfin::end ();
    }



    void InstantaneousControl::solve (const std::string& problemName, const bool& forceRelinking)
    {
        solve (problemName, forceRelinking, true);
    }



    void InstantaneousControl::solve (const char* problemName, const bool& forceRelinking)
    {
        solve (std::string (problemName), forceRelinking);
    }
    

    void InstantaneousControl::setMeshManager (const MeshManager<dolfin::ALE> & meshManager)
      { meshManager_ = & meshManager; }

    dcp::MeshManager<> InstantaneousControl::meshManager () const
      { return * meshManager_; }

    /******************* PROTECTED METHODS *******************/
    void InstantaneousControl::solve (const std::string& problemName, const bool& forceRelinking, const bool advanceTimeDependentProblem)
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
        
        auto tdProblemIterator = timeDependentStoredProblems_.find (problemName);
        if (tdProblemIterator != timeDependentStoredProblems_.end())
          if (advanceTimeDependentProblem)
            problem.solve ("step");
          else
            problem.solve ("noAdvance_step");
        else
          problem.solve ();
        
        dolfin::end ();
    }

    double InstantaneousControl::chooseIncrementStep (const double alpha, const double alphaMin, dolfin::Array<double> & controlValues, const dolfin::Array<double> & controlPoint, const std::string & controlName, const double & projectedAdjoint, const bool& forceRelinking)
    {
        if (alpha<=alphaMin)
          return alphaMin;

        control_->eval (controlValues, controlPoint);
        double rescaleValue (rescale_.value());
        /*control_.reset (new dolfin::Constant (0, 
            controlValues[1]
            - alpha
            * static_cast<dcp::TimeDependentProblem&> ((*this)["primal"]).dt()
            * ( controlValues[1] * rescaleValue * problemData.lx
//                * static_cast<dcp::TimeDependentProblem&> ((*this)["primal"]).time()
//                * static_cast<dcp::TimeDependentProblem&> ((*this)["primal"]).time()
                + projectedAdjoint )
            ) );  */
        control_.reset (new dolfin::Constant (0, 
            controlValues[1]
            + alpha
            * static_cast<dcp::TimeDependentProblem&> ((*this)["primal"]).dt()
            * ( projectedAdjoint
                - rescaleValue * problemData.lx * controlValues[1] ) ) );

        (*this)["primal"].setCoefficient ("linear_form", control_, controlName);
        solve("primal", forceRelinking, false);

        objectiveFunctional_->set_coefficient ("u", dolfin::reference_to_no_delete_pointer ((*this)["primal"].solution()[0]));
        objectiveFunctional_->set_coefficient (controlName, control_);
        //rescale_.setTime (static_cast<dcp::TimeDependentProblem&> ((*this)["primal"]).time());
        double objFunVal (dolfin::assemble (* objectiveFunctional_));
std::cerr << "VALS " << objFunVal << "  " << objFunValues_.back () << "   alpha " << alpha << std::endl;

        if (objFunVal < objFunValues_.back ())
        {
            objFunValues_.push_back (dolfin::assemble (* objectiveFunctional_));
            return alpha;
        }
        else
            return this->chooseIncrementStep (0.5 * alpha, alphaMin, controlValues, controlPoint, controlName, projectedAdjoint, forceRelinking);
    }

    void InstantaneousControl::linkProblemToPreviousSolution (const PreviousSolutionLink& link)
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
/*std::cerr << __LINE__ << " ; ";
std::cerr << targetProblemSolutionsVector.size() << " ; ";
std::cerr << std::boolalpha << (targetProblemSolutionsVector.rend() == (targetProblemSolutionsVector.rbegin() + nStepsBack)) << " ; ";
std::cerr << (targetProblemSolutionsVector.rend() - (targetProblemSolutionsVector.rbegin() + nStepsBack)) << " ; ";
std::cerr << targetFunction.id() << " ; ";
std::cerr << std::get<0> (link.second) << " ; ";
std::cerr << targetFunction.str(true) << " ; ";
std::cerr << component << std::endl;*/
            problem.setCoefficient (std::get<2> (link.first), 
                                    dolfin::reference_to_no_delete_pointer (targetFunction [component]),
                                    std::get<1> (link.first));

        }
        
        dolfin::end ();
    }
}
