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

//#include <dcp/differential_problems/MovingTimeDependentProblem.h>
#include "MovingTimeDependentProblem.h"
#include "computeFreeSurfaceStress_onlyTP.h"
#include "computeFreeSurfaceStress_noTP.h"
#include <dcp/differential_problems/differential_problems.h>
#include <dolfin/log/dolfin_log.h>
#include <dolfin/geometry/Point.h>
#include <regex>
#include "UflToNewton.h"
#include "myNavierstokesTimeCurv.h"

dolfin::Point extractPoint (dolfin::MeshFunction<std::size_t>& mf, size_t value)
{
    std::vector<std::size_t> vec;
    vec.assign (mf.values(), mf.values()+mf.size());
    std::vector<std::size_t>::iterator itP1 = find (vec.begin(),vec.end(),value);
    std::size_t idxP1 = std::distance (vec.begin(), itP1);
    return (mf.mesh()->geometry().point(idxP1));
}

namespace Ivan
{
    /******************* CONSTRUCTORS *******************/
    MovingTimeDependentProblem::MovingTimeDependentProblem 
        (const std::shared_ptr<dcp::AbstractProblem> timeSteppingProblem,
         const double& startTime,
         const double& dt,
         const double& endTime,
         std::initializer_list<std::string> dtCoefficientTypes,
         std::initializer_list<std::string> previousSolutionCoefficientTypes,
//current//         std::initializer_list<std::string> currentSolutionCoefficientTypes,
         const unsigned int& nTimeSchemeSteps)
        : 
            AbstractProblem (timeSteppingProblem->functionSpace ()),
            timeSteppingProblem_ (timeSteppingProblem),
            timeDependentCoefficients_ (),
            timeDependentDirichletBCs_ (),
            t_ (startTime),
            startTime_ (startTime),
            dt_ (dt),
            endTime_ (endTime),
            nTimeSchemeSteps_ (nTimeSchemeSteps),
            timeDependentDirichletBCsCounter_ (0),
            meshManager_ (std::shared_ptr<dolfin::ALE>(new dolfin::ALE()),timeSteppingProblem->functionSpace())
    { 
        dolfin::begin (dolfin::DBG, "Building MovingTimeDependentProblem...");
        
        // we need a solution for every step in the time scheme! 
        // And the time should be t0, t0+dt, t0+2dt ...
        dolfin::begin (dolfin::DBG, "Creating initial solutions...");
        unsigned int stepNumber;
        for (stepNumber = 0; stepNumber < nTimeSchemeSteps_; ++stepNumber)
        {
            solution_.emplace_back (std::make_pair (t_ + stepNumber * dt, 
                                                    dolfin::Function (timeSteppingProblem_->functionSpace ())));
        }
        dolfin::log (dolfin::DBG, "Created %d initial solutions", stepNumber + 1);
        
        // set the correct value for t_, since it was not incremented during the previous loop.
        // Use stepNumber - 1 since t_ will be incremented at the *beginning* of the time loop, not at the end
        t_ += (stepNumber - 1) * dt;
            
        dolfin::end ();
        
        dolfin::log (dolfin::DBG, "Setting up parameters...");
        parameters.add ("problem_type", "time_dependent");
        parameters.add ("dt_name", "dt");
        parameters.add ("previous_solution_name", "u_old");
        parameters.add ("current_solution_name", "trial");
        parameters.add ("store_interval", 1);
        parameters.add ("plot_interval", 0);
//        parameters.add ("time_stepping_solution_component", -1);
        parameters.add ("time_stepping_solution_component", 0);
          // TODO rimettere default a -1 e capire come passarlo da fuori
        parameters.add ("time_stepping_current_solution_component", -1);
        parameters.add ("pause", false);
        parameters.add ("previous_solution_is_set_externally", false);
        
        dolfin::Parameters dtCoefficientTypesParameter ("dt_coefficient_types");
        for (auto& i : dtCoefficientTypes)
        {
            dtCoefficientTypesParameter.add<std::string> (i);
        }
        parameters.add (dtCoefficientTypesParameter);
        
        dolfin::Parameters previousSolutionCoefficientTypesParameter ("previous_solution_coefficient_types");
        for (auto& i : previousSolutionCoefficientTypes)
        {
            previousSolutionCoefficientTypesParameter.add<std::string> (i);
        }
        parameters.add (previousSolutionCoefficientTypesParameter);
        
//current//        dolfin::Parameters currentSolutionCoefficientTypesParameter ("current_solution_coefficient_types");
//current//        for (auto& i : currentSolutionCoefficientTypes)
//current//        {
//current//            currentSolutionCoefficientTypesParameter.add<std::string> (i);
//current//        }
//current//        parameters.add (currentSolutionCoefficientTypesParameter);

        parameters ["plot_title"] = "";

        dolfin::end ();
        
        dolfin::log (dolfin::DBG, "MovingTimeDependentProblem object created");
    }



    /******************* GETTERS *******************/
    std::shared_ptr<const dolfin::Mesh> MovingTimeDependentProblem::mesh () const
    {
//        return timeSteppingProblem_ -> mesh ();
        return dolfin::reference_to_no_delete_pointer(meshManager_.mesh ());
    }



    std::shared_ptr<dolfin::FunctionSpace> MovingTimeDependentProblem::functionSpace () const
    {
        return timeSteppingProblem_ -> functionSpace ();
    }



    const dolfin::DirichletBC& MovingTimeDependentProblem::dirichletBC (const std::string& bcName) const
    {
        return timeSteppingProblem_ -> dirichletBC (bcName);
    }

    

    const std::map<std::string, dolfin::DirichletBC>& MovingTimeDependentProblem::dirichletBCs () const
    {
        return timeSteppingProblem_ -> dirichletBCs ();
    }
    


    const dolfin::Function& MovingTimeDependentProblem::solution () const
    {
        return solution_.back ().second;
    }



    const std::vector <std::pair <double, dolfin::Function> >& MovingTimeDependentProblem::solutions () const
    {
        return solution_;
    }
    


    const double& MovingTimeDependentProblem::time () const
    {
        return t_;
    }
    


    double& MovingTimeDependentProblem::startTime ()
    {
        return startTime_;
    }
    


    double& MovingTimeDependentProblem::dt ()
    {
        return dt_;
    }



    double& MovingTimeDependentProblem::endTime ()
    {
        return endTime_;
    }
    


    dcp::AbstractProblem& MovingTimeDependentProblem::timeSteppingProblem ()
    {
        return *timeSteppingProblem_;
    }



    /******************* SETTERS *******************/
    void MovingTimeDependentProblem::setInitialSolution (const dolfin::Function& initialSolution, 
                                                   const unsigned int& stepNumber)
    {
        if (stepNumber > nTimeSchemeSteps_)
        {
            dolfin::warning ("initial solution not set. Requested time step number is greater than problem's scheme's number of time steps");
            return;
        }
        
        auto solutionsIterator = solution_.rbegin ();
        
        // the solution we want to set is identified by stepNumber. If it is equal to 1, we must set the last element
        // in solution_, if it is equal to 2 the last but one element and so on. That's why we subtract 1 from 
        // stepNumber in the next instruction
        (solutionsIterator + (stepNumber - 1)) -> second = initialSolution;
    }

    

    void MovingTimeDependentProblem::setInitialSolution (const dolfin::Expression& initialSolution,
                                                   const unsigned int& stepNumber)
    {
        if (stepNumber > nTimeSchemeSteps_)
        {
            dolfin::warning ("initial solution not set. Requested time step number is greater than problem's scheme's number of time steps");
            return;
        }
        
        auto solutionsIterator = solution_.rbegin ();
        
        // the solution we want to set is identified by stepNumber. If it is equal to 1, we must set the last element
        // in solution_, if it is equal to 2 the last but one element and so on. That's why we subtract 1 from 
        // stepNumber in the next instruction
        (solutionsIterator + (stepNumber - 1)) -> second = initialSolution;
    }
   


    void MovingTimeDependentProblem::setCoefficient (const std::string& coefficientType, 
                                               const std::shared_ptr<const dolfin::GenericFunction> coefficientValue,
                                               const std::string& coefficientName)
    {
        timeSteppingProblem_->setCoefficient (coefficientType, coefficientValue, coefficientName);
    }



    void MovingTimeDependentProblem::setCoefficient (const std::string& coefficientType,
                                               const std::shared_ptr<const dolfin::GenericFunction> coefficientValue,
                                               const std::size_t& coefficientNumber)
    {
        timeSteppingProblem_->setCoefficient (coefficientType, coefficientValue, coefficientNumber);
    }



    void MovingTimeDependentProblem::
    setIntegrationSubdomain (const std::string& formType,
                              std::shared_ptr<const dolfin::MeshFunction<std::size_t>> meshFunction,
                              const dcp::SubdomainType& subdomainType)
    {
        timeSteppingProblem_->setIntegrationSubdomain (formType, meshFunction, subdomainType);
    }



    bool MovingTimeDependentProblem::addDirichletBC (const dolfin::GenericFunction& condition, 
                                               const dolfin::SubDomain& boundary,
                                               std::string bcName)
    {
        return timeSteppingProblem_->addDirichletBC (condition, boundary, bcName);
    }
    


    bool MovingTimeDependentProblem::addDirichletBC (const dolfin::GenericFunction& condition, 
                                               const dolfin::SubDomain& boundary, 
                                               const std::size_t& component,
                                               std::string bcName)
    {
        return timeSteppingProblem_->addDirichletBC (condition, boundary, component, bcName);
    }
    


    bool MovingTimeDependentProblem::addDirichletBC (std::shared_ptr<const dolfin::GenericFunction> condition, 
                                               std::shared_ptr<const dolfin::SubDomain> boundary,
                                               std::string bcName)
    {
        return timeSteppingProblem_->addDirichletBC (condition, boundary, bcName);
    }
    


    bool MovingTimeDependentProblem::addDirichletBC (std::shared_ptr<const dolfin::GenericFunction> condition, 
                                               std::shared_ptr<const dolfin::SubDomain> boundary, 
                                               const std::size_t& component,
                                               std::string bcName)
    {
        return timeSteppingProblem_->addDirichletBC (condition, boundary, component, bcName);
    }
    


    bool MovingTimeDependentProblem::addDirichletBC (const dolfin::DirichletBC& dirichletCondition, 
                                               std::string bcName)
    {
        return timeSteppingProblem_->addDirichletBC (dirichletCondition, bcName);
    }



    bool MovingTimeDependentProblem::addDirichletBC (dolfin::DirichletBC&& dirichletCondition, 
                                               std::string bcName)
    {
        return timeSteppingProblem_->addDirichletBC (dirichletCondition, bcName);
    }



    bool MovingTimeDependentProblem::removeDirichletBC (const std::string& bcName)
    {
        return timeSteppingProblem_->removeDirichletBC (bcName);
    }

    
        
    bool MovingTimeDependentProblem::addTimeDependentDirichletBC (const dcp::TimeDependentExpression& condition, 
                                                            const dcp::Subdomain& boundary,
                                                            std::string bcName)
    {
        if (bcName.empty ())
        {
            bcName = "time_dependent_dirichlet_condition_" + std::to_string (timeDependentDirichletBCsCounter_);
            timeDependentDirichletBCsCounter_++;
        }

        dolfin::begin (dolfin::DBG, 
                       "Adding dirichlet boundary condition to time dependent boundary conditions map with name \"%s\"...",
                       bcName.c_str ());

        auto result = timeDependentDirichletBCs_.emplace 
            (bcName, 
             std::make_tuple (std::shared_ptr<dcp::TimeDependentExpression> (condition.clone ()), 
                              std::shared_ptr<dcp::Subdomain> (boundary.clone ()), 
                              -1)
             );

        bool bcWasAdded = result.second;

        if (bcWasAdded == false)
        {
            dolfin::warning ("DirichletBC object not inserted because key \"%s\" already in map",
                             bcName.c_str ());
        }
        else
        {
            dolfin::begin (dolfin::DBG, 
                           "Adding dirichlet boundary condition to time stepping problem boundary conditions map with name \"%s\"...",
                           bcName.c_str ());
            bcWasAdded = bcWasAdded && addDirichletBC (std::get<0> ((result.first)->second),
                                                       std::get<1> ((result.first)->second),
                                                       bcName);
            dolfin::end ();
        }
        
        dolfin::end ();

        return bcWasAdded;
    }



    bool MovingTimeDependentProblem::addTimeDependentDirichletBC (const dcp::TimeDependentExpression& condition, 
                                                            const dcp::Subdomain& boundary,
                                                            const std::size_t& component,
                                                            std::string bcName)
    {
        if (bcName.empty ())
        {
            bcName = "time_dependent_dirichlet_condition_" + std::to_string (timeDependentDirichletBCsCounter_);
            timeDependentDirichletBCsCounter_++;
        }

        dolfin::begin (dolfin::DBG, 
                       "Adding dirichlet boundary condition to time dependent boundary conditions map with name \"%s\"...",
                       bcName.c_str ());

        auto result = timeDependentDirichletBCs_.emplace 
            (bcName, 
             std::make_tuple (std::shared_ptr<dcp::TimeDependentExpression> (condition.clone ()), 
                              std::shared_ptr<dcp::Subdomain> (boundary.clone ()), 
                              component)
             );

        bool bcWasAdded = result.second;

        if (bcWasAdded == false)
        {
            dolfin::warning ("DirichletBC object not inserted because key \"%s\" already in map",
                             bcName.c_str ());
        }
        else
        {
            dolfin::begin (dolfin::DBG, 
                           "Adding dirichlet boundary condition to time stepping problem boundary conditions map with name \"%s\"...",
                           bcName.c_str ());
            bcWasAdded = bcWasAdded && addDirichletBC (std::get<0> ((result.first)->second),
                                                       std::get<1> ((result.first)->second),
                                                       component,
                                                       bcName);
            dolfin::end ();
        }

        dolfin::end ();

        return result.second;
    }
    
    

    bool MovingTimeDependentProblem::removeTimeDependentDirichletBC (const std::string& bcName)
    {
        dolfin::begin (dolfin::DBG, 
                     "Removing dirichlet boundary condition \"%s\" from time dependent boundary conditions map...", 
                     bcName.c_str ());
        std::size_t nErasedElements = timeDependentDirichletBCs_.erase (bcName);

        if (nErasedElements == 0)
        {
            dolfin::warning ("Dirichlet boundary condition \"%s\" not found in map", 
                             bcName.c_str ());
        }
        
        dolfin::begin (dolfin::DBG, 
                       "Removing dirichlet boundary condition \"%s\" from time stepping problem boundary conditions map...", 
                       bcName.c_str ());
        // call removeDirichletBC to remove the dirichlet bc also from timeSteppingProblem_
        // We add it to the existing value of nErasedElements so that we know wether both removal operations went fine
        nErasedElements += removeDirichletBC (bcName);
        dolfin::end ();

        dolfin::end ();
        
        // remember: both removal were ok if nErasedElements was equal to 1 both times!
        return nErasedElements == 2? true : false;
    }



    bool MovingTimeDependentProblem::addTimeDependentCoefficient (const std::string& coefficientName, 
                                                            const std::string& coefficientType,
                                                            std::shared_ptr<dcp::TimeDependentExpression> expression)
    {
        dolfin::begin (dolfin::DBG, 
                       "Inserting time dependent coefficient in map with name \"%s\" and type \"%s\"...",
                       coefficientName.c_str (),
                       coefficientType.c_str ());
       
        auto coefficientID = std::make_pair (coefficientName, coefficientType);
        auto result = timeDependentCoefficients_.insert (std::make_pair (coefficientID, expression));
        
        if (result.second == false)
        {
            dolfin::warning ("Time dependent coefficient not inserted because key is already present in map");
        }
        
        dolfin::end ();
        
        return result.second;
    }



    bool MovingTimeDependentProblem::removeTimeDependentCoefficient (const std::string& coefficientName,
                                                               const std::string& coefficientType)
    {
        dolfin::begin (dolfin::DBG, 
                       "Removing time dependent coefficient with name \"%s\" and type \"%s\" from map...",
                       coefficientName.c_str (),
                       coefficientType.c_str ());
       
        auto coefficientID = std::make_pair (coefficientName, coefficientType);
        
        std::size_t nErasedElements = timeDependentCoefficients_.erase (coefficientID);
        
        if (nErasedElements == 0)
        {
            dolfin::warning ("Time dependent coefficient not removed: key was not found in map");
        }
        
        dolfin::end ();
        
        return nErasedElements == 1? true : false;
    }
            


    void MovingTimeDependentProblem::update ()
    {
        timeSteppingProblem_->update ();
    }



    /******************* METHODS *******************/
    bool MovingTimeDependentProblem::isFinished ()
    {
        return (dt_ > 0) ? (t_ >= endTime_ - DOLFIN_EPS) : (t_ <= endTime_ + DOLFIN_EPS);
    }
    


    void MovingTimeDependentProblem::clear ()
    {
        // reset t_
        t_ = startTime_;
        
        // clear solutions vector
        solution_.clear ();
        solution_.emplace_back (std::make_pair (t_, dolfin::Function (timeSteppingProblem_->functionSpace ())));
        
        // reset time dependent Dirichlet BCs
        for (auto bcIterator = timeDependentDirichletBCs_.begin (); 
             bcIterator != timeDependentDirichletBCs_.end (); 
             bcIterator++)
        {
            resetTimeDependentDirichletBC (bcIterator);
        }
    }
    


    void MovingTimeDependentProblem::step ()
    {
        solve ("step");
    }
    

    
    void MovingTimeDependentProblem::solve (const std::string& type) 
    {
        // parse solve type
        if (type != "default" && type != "step" && type != "clear_default" && type != "clear_step")
        {
            dolfin::dolfin_error ("dcp: MovingTimeDependentProblem.h", 
                                  "solve",
                                  "Unknown solve type \"%s\" requested",
                                  type.c_str ());
        }
        
        // check if we are already at the end of time loop
        if (isFinished ())
        {
            printFinishedWarning ();
            return;
        }

        dolfin::log (dolfin::DBG, "Solve type: %s", type.c_str ());
            
        if (std::regex_match (type, std::regex (".*clear.*")))
        {
            clear ();
        }
        
        // ---- Problem settings ---- //
        dolfin::begin (dolfin::DBG, "Setting up time dependent problem...");
        
        // flag to check if the solver must perform just one time step
        bool oneStepRequested = std::regex_match (type, std::regex (".*step.*"));
        
        // get parameters' values
        dolfin::Constant dt (dt_);
        std::string dtName = parameters ["dt_name"];
        std::vector<std::string> dtCoefficientTypes;
        parameters ("dt_coefficient_types").get_parameter_keys (dtCoefficientTypes);
        int storeInterval = parameters ["store_interval"];
        int plotInterval = parameters ["plot_interval"];
        int plotComponent = parameters ["plot_component"];
        bool pause = parameters ["pause"];
        
        for (auto& i : dtCoefficientTypes)
        {
            timeSteppingProblem_->setCoefficient (i, dolfin::reference_to_no_delete_pointer (dt), dtName);
        } 
        
        dolfin::end ();
        
        
        // ---- Problem solution ---- //
        dolfin::begin (dolfin::DBG, "Solving time dependent problem...");
        
        
        // function used to step through the time loop
        dolfin::Function tmpSolution = solution_.back ().second;
        
        std::shared_ptr<dolfin::VTKPlotter> plotter;
        
//lasciarlo qui?? l'importante è che il marking avvenga solo una volta (quindi prima del ciclo temporale)
// è per calcolare la tensione superficiale
    dolfin::FacetFunction<std::size_t> meshFacets (meshManager_.mesh());
    meshFacets.set_all (0);
    TopBoundary freeSurface;
    freeSurface.mark (meshFacets, 1);
    LateralWall wallBoundary;
    wallBoundary.mark (meshFacets, 2);
Ivan::TriplePointLeft triplePointLeft;
triplePointLeft.mark(meshFacets,4);
Ivan::TriplePointRight triplePointRight;
triplePointRight.mark(meshFacets,5);
dolfin::plot(meshFacets,"MovingTimeDepPb meshFacets");//dolfin::interactive();
/*dolfin::VertexFunction<std::size_t> meshPoints (meshManager_.mesh());
meshPoints.set_all (0);
triplePointLeft.mark(meshPoints,4);
triplePointRight.mark(meshPoints,5);
meshPoints.set_value(220,10);
meshPoints.set_value(230,20);
dolfin::Point tpLeft (extractPoint(meshPoints,10));
dolfin::Point tpRight (extractPoint(meshPoints,20));
dolfin::plot(meshPoints, "meshPoints"); dolfin::interactive();*/
/*this->setIntegrationSubdomain ("residual_form",
                                            dolfin::reference_to_no_delete_pointer (meshFacets),
                                            dcp::SubdomainType::BOUNDARY_FACETS);*/
/* delta// computeFreeSurfaceStress::Form_L::CoefficientSpace_deltaDirac deltaFS (this->mesh());
//dolfin::PointSource deltaDirac (deltaFS,tpLeft);
Ivan::DeltaDirac deltaExpr;
dolfin::Function deltaFun (deltaFS);
deltaFun = deltaExpr;
dolfin::plot(deltaFun,"la delta fun"); dolfin::interactive();
dolfin::plot(deltaExpr,meshManager_.mesh(),"la delta expr");
delta// */
											
        // start time loop. The loop variable timeStepFlag is true only if the solve type requested is "step" and
        // the first step has already been peformed
        int timeStep = 0;
        bool timeStepFlag = false;
dolfin::File solutionFile("/u/laureandi/ifumagalli/dcp_test_output/insideMovingTimeDependentProblem.pvd");
solutionFile << std::pair<dolfin::Function*,double>(&tmpSolution,t_);
dolfin::File surfaceStressFile("/u/laureandi/ifumagalli/dcp_test_output/surfaceTension.pvd");
dolfin::File surfaceStressFile_onlyTP("/u/laureandi/ifumagalli/dcp_test_output/surfaceTension_onlyTP.pvd");
std::ofstream surfaceStressFile_onlyTP_csv; surfaceStressFile_onlyTP_csv.open("/u/laureandi/ifumagalli/dcp_test_output/surfaceTension_onlyTP.csv",std::ios::out|std::ios::app);
dolfin::File surfaceStressFile_noTP("/u/laureandi/ifumagalli/dcp_test_output/surfaceTension_noTP.pvd");
dolfin::File provaFile("/u/laureandi/ifumagalli/dcp_test_output/provaFile.pvd");
        while (!isFinished () && timeStepFlag == false)
        {
            timeStep++;

            meshManager_.moveMesh(tmpSolution[0],"normal",dt);
            dolfin::Function& w(meshManager_.displacement());
//dolfin::plot(w,"w in MovingTimeDependentProblem");
//            dolfin::Constant w(0.0,0.0001);
//dolfin::plot(w,meshManager_.mesh(),"w in MovingTimeDependentProblem");

//dolfin::interactive();
//            timeSteppingProblem_->setIntegrationSubdomain ("residual_form", meshFunction, dcp::SubdomainType::BOUNDARY_FACETS);
            timeSteppingProblem_->setCoefficient ("residual_form", dolfin::reference_to_no_delete_pointer (w), "w");
            timeSteppingProblem_->setCoefficient ("jacobian_form", dolfin::reference_to_no_delete_pointer (w), "w");
            
            if (oneStepRequested == false)
            {
                dolfin::begin (dolfin::INFO, "===== Timestep %d =====", timeStep);
            }
            
            advanceTime ();
            
            setTimeDependentCoefficients ();
            
            setTimeDependentDirichletBCs ();
            
            dolfin::begin (dolfin::INFO, "Solving time stepping problem...");
provaFile << timeSteppingProblem_->solution ();
            timeSteppingProblem_->solve ();
            dolfin::end ();
            
            tmpSolution = timeSteppingProblem_->solution ();
provaFile << tmpSolution;
dolfin::Function newTimeStep (* tmpSolution.function_space());
newTimeStep = dolfin::Constant (-0.001,0,0);
provaFile << newTimeStep;
solutionFile << std::pair<dolfin::Function*,double>(&tmpSolution,t_);
            
// Assembling and storing the surface tension term
//delta// dolfin::plot(deltaFun,"la delta fun"); 
computeFreeSurfaceStress_noTP::FunctionSpace stressSpace_noTP(meshManager_.mesh());
computeFreeSurfaceStress_onlyTP::FunctionSpace stressSpace(meshManager_.mesh());
dolfin::Constant stressGamma (7.3e-5);
dolfin::Constant stressDt (0.05);
// CI FOSSE residualForm() COME METODO PER ESTRARRE I COEFFICIENTI...  stressForm.set_coefficient ("gamma",timeSteppingProblem_->residualForm().coefficient("gamma"));
//delta// computeFreeSurfaceStress::LinearForm stressForm(stressSpace,stressDt,stressGamma);//,deltaFun);
dolfin::Constant beta (1.0e-5);
dolfin::Constant cosThetaS (cos(3.14159265/3.0));
dolfin::Constant t_partialOmega(0,1);
computeFreeSurfaceStress_noTP::LinearForm stressForm_noTP(stressSpace_noTP,tmpSolution,stressDt,stressGamma,beta);//,deltaFun);
computeFreeSurfaceStress_onlyTP::LinearForm stressForm_onlyTP(stressSpace,stressDt,stressGamma,cosThetaS,t_partialOmega);//,deltaFun);
//delta// stressForm.set_exterior_facet_domains (dolfin::reference_to_no_delete_pointer(meshFacets));
stressForm_noTP.set_exterior_facet_domains (dolfin::reference_to_no_delete_pointer(meshFacets));
stressForm_onlyTP.set_exterior_facet_domains (dolfin::reference_to_no_delete_pointer(meshFacets));
//stressForm.set_vertex_domains (dolfin::reference_to_no_delete_pointer(meshPoints));
dolfin::Vector b;
dolfin::Function bFun(stressSpace);
/* delta// assemble(b, stressForm);            
* bFun.vector() = b;
//dolfin::plot(bFun,"bFun"); dolfin::interactive();
surfaceStressFile << std::pair<dolfin::Function*,double>(&bFun,t_);
delta// */
assemble(b, stressForm_onlyTP);            
* bFun.vector() = b;
surfaceStressFile_onlyTP << std::pair<dolfin::Function*,double>(&bFun,t_);
for (std::size_t i=0; i!=b.size(); ++i) surfaceStressFile_onlyTP_csv << b[i] << ","; surfaceStressFile_onlyTP_csv << std::endl;
dolfin::Function bFun_noTP(stressSpace_noTP);
assemble(* bFun_noTP.vector(), stressForm_noTP);
surfaceStressFile_noTP << std::pair<dolfin::Function*,double>(&bFun_noTP,t_);
std::cerr << std::endl;

            // save solution in solution_ according to time step and store interval.
            storeSolution (tmpSolution, timeStep, storeInterval);
            
            // plot solution according to time step and plot interval
 //           plotSolution (tmpSolution, timeStep, plotInterval, plotComponent, pause);
            plotSolution (tmpSolution, timeStep, 1, 0, false);//dolfin::interactive();

            if (oneStepRequested == true)
            {
                timeStepFlag = true;
            }
             
            // dolfin::end matching dolfin::begin inside if statement
            if (oneStepRequested == false)
            {
                dolfin::end ();
            }
        
        }
        
        // At this point, we just need to make sure that the solution on the last iteration was saved even though
        // timeStep % storeInterval != 0 (but it must not be saved twice!)
        storeLastStepSolution (tmpSolution, timeStep, storeInterval);

surfaceStressFile_onlyTP_csv.close();
        dolfin::end ();
    }
    
    

    void MovingTimeDependentProblem::plotSolution ()
    {
        bool pause = parameters ["pause"];
        int plotComponent = parameters ["plot_component"];
        std::string plotTitle = parameters ["plot_title"];
        
        // if plotTitle is not empty, we need to prepend ", " so that the plot title is readable
        if (!plotTitle.empty ())
        {
            plotTitle = ", " + plotTitle;
        }
        
        // auxiliary variable, to enhance readability
        std::shared_ptr<dolfin::Function> functionToPlot;
        
        dolfin::begin (dolfin::DBG, "Plotting...");
        
        for (auto& timeSolutionPair : solution_)
        {
            double time = timeSolutionPair.first;
            // get right function to plot
            if (plotComponent == -1)
            {
                functionToPlot = dolfin::reference_to_no_delete_pointer (timeSolutionPair.second);
                dolfin::log (dolfin::DBG, "Plotting time stepping problem solution, all components...");
            }
            else
            {
                functionToPlot = dolfin::reference_to_no_delete_pointer (timeSolutionPair.second [plotComponent]);
                dolfin::log (dolfin::DBG, "Plotting time stepping problem solution, component %d...", plotComponent);
            }

            // actual plotting
            if (solutionPlotter_ == nullptr)
            {
                dolfin::log (dolfin::DBG, "Plotting in new dolfin::VTKPlotter object...");
                solutionPlotter_ = dolfin::plot (functionToPlot, "Time = " + std::to_string (time) + plotTitle);
            }
            else if (! solutionPlotter_ -> is_compatible (functionToPlot))
            {
                dolfin::log (dolfin::DBG, "Existing plotter is not compatible with object to be plotted.");
                dolfin::log (dolfin::DBG, "Creating new dolfin::VTKPlotter object...");
                solutionPlotter_ = dolfin::plot (functionToPlot, "Time = " + std::to_string (time) + plotTitle);
            }
            else 
            {
                solutionPlotter_ -> parameters ["title"] = std::string ("Time = " + std::to_string (time) + plotTitle);
                solutionPlotter_ -> plot (functionToPlot);
            }

            if (pause)
            {
                dolfin::interactive ();
            }
        }
        
        dolfin::end ();
    }



    Ivan::MovingTimeDependentProblem* MovingTimeDependentProblem::clone () const
    {
        dolfin::begin (dolfin::DBG, "Cloning object...");
        
        std::string cloneMethod = parameters ["clone_method"];
        
        dolfin::log (dolfin::DBG, "Clone method: %s", cloneMethod.c_str ());
        dolfin::log (dolfin::DBG, "Creating new object of type MovingTimeDependentProblem...");
        
        // create new object
        Ivan::MovingTimeDependentProblem* clonedProblem = nullptr;
        if (cloneMethod == "shallow_clone")
        {
            // note that we pass an empty initializer_list to the constructor as dtCoefficientTypes and 
            // previousSolutionCoefficientTypes, because they will be copied when the parameters are copied anyway
            clonedProblem = 
                new Ivan::MovingTimeDependentProblem (this->timeSteppingProblem_, startTime_, dt_, endTime_, {}, {}, {});
            clonedProblem->timeDependentDirichletBCs_ = this->timeDependentDirichletBCs_;
        }
        else if (cloneMethod == "deep_clone")
        {
            // note that we pass an empty initializer_list to the constructor as dtCoefficientTypes and 
            // previousSolutionCoefficientTypes, because they will be copied when the parameters are copied anyway
            clonedProblem = 
                new Ivan::MovingTimeDependentProblem (this->timeSteppingProblem_, startTime_, dt_, endTime_, {}, {}, {});
        }
        else
        {
            dolfin::dolfin_error ("dcp: MovingTimeDependentProblem.cpp",
                                  "clone",
                                  "Cannot clone time dependent differential problem. Unknown clone method: \"%s\"",
                                  cloneMethod.c_str ());
            for (auto& bc : this->timeDependentDirichletBCs_)
            {
                clonedProblem->addTimeDependentDirichletBC (*(std::get<0> (bc.second)),
                                                            *(std::get<1> (bc.second)),
                                                            std::get<2> (bc.second),
                                                            bc.first);
            }
        }
        
        //copy dirichlet boundary conditions
        dolfin::log (dolfin::DBG, "Copying Dirichlet boundary conditions...");
        for (auto& bc : this->dirichletBCs_)
        {
            clonedProblem->addDirichletBC (bc.second, bc.first);
        }

        // clear parameters set of newly created object so that it can be populated by the parameters of the object
        // being created. 
        dolfin::log (dolfin::DBG, "Copying parameters to new object...");
        clonedProblem->parameters.clear ();
        clonedProblem->parameters = this->parameters;
        
        // copy solution
        dolfin::log (dolfin::DBG, "Copying solution...");
        clonedProblem->solution_ = this->solution_;
        
        dolfin::end ();
        
        return clonedProblem;
    }
    


    /******************* PROTECTED METHODS *******************/
    void MovingTimeDependentProblem::advanceTime ()
    {
        t_ += dt_;
        
        dolfin::log (dolfin::INFO, "TIME = %f s", t_);

        std::string previousSolutionName = parameters ["previous_solution_name"];
        std::vector<std::string> previousSolutionCoefficientTypes;
        parameters ("previous_solution_coefficient_types").get_parameter_keys (previousSolutionCoefficientTypes);
        int timeSteppingSolutionComponent = parameters ["time_stepping_solution_component"];
        
        bool previousSolutionIsSetExternally = parameters ["previous_solution_is_set_externally"];
        if (previousSolutionIsSetExternally == false)
        {
            dolfin::begin (dolfin::DBG, "Setting previous solution coefficients...");
            for (unsigned int step = 1; step <= nTimeSchemeSteps_; ++step)
            {
                // if step > 1, append step number to the name of the coefficient to set 
                if (step > 1)
                {
                    previousSolutionName = previousSolutionName + "_" + std::to_string (step);
                }
                
                // get the index of the function in solution_ to use to set the coefficient
                unsigned int index = solution_.size () - step;
                for (auto& previousSolutionCoefficientType : previousSolutionCoefficientTypes)
                {
std::cerr << previousSolutionName << " : COMPONENTE " << timeSteppingSolutionComponent << std::endl;
                    if (timeSteppingSolutionComponent >= 0)
                    {
                        timeSteppingProblem_ -> setCoefficient 
                            (previousSolutionCoefficientType, 
                             dolfin::reference_to_no_delete_pointer ((solution_ [index].second) [timeSteppingSolutionComponent]), 
                             previousSolutionName);
                    }
                    else
                    {
                        timeSteppingProblem_ -> setCoefficient 
                            (previousSolutionCoefficientType, 
                             dolfin::reference_to_no_delete_pointer (solution_ [index].second),
                             previousSolutionName);
                    }
                } 
            }
            dolfin::end ();
        }
        else
        {
            dolfin::log (dolfin::DBG, "Skipping previous solution setting loop since it is set externally.");
        }

/*current//        std::string currentSolutionName = parameters ["current_solution_name"];
        std::vector<std::string> currentSolutionCoefficientTypes;
        parameters ("current_solution_coefficient_types").get_parameter_keys (currentSolutionCoefficientTypes);
        int timeSteppingCurrentSolutionComponent = parameters ["time_stepping_current_solution_component"];

        dolfin::begin (dolfin::DBG, "Setting current solution as coefficient (nonlinear solver)...");
          // TODO bisogna vedere meglio dove mettere questa cosa...
        for (auto& currentSolutionCoefficientType : currentSolutionCoefficientTypes)
        {
std::cerr << currentSolutionName << " : COMPONENTE " << timeSteppingCurrentSolutionComponent << std::endl;
            if (timeSteppingCurrentSolutionComponent >= 0)
            {
                timeSteppingProblem_ -> setCoefficient 
                    (currentSolutionCoefficientType, 
                     dolfin::reference_to_no_delete_pointer ((solution_.back().second) [timeSteppingCurrentSolutionComponent]), 
                     currentSolutionName);
            }
            else
            {
                timeSteppingProblem_ -> setCoefficient 
                    (currentSolutionCoefficientType, 
                     dolfin::reference_to_no_delete_pointer (solution_.back().second),
                     currentSolutionName);
            }
        }
        dolfin::end ();
//current*/
    }
            
    

    void MovingTimeDependentProblem::setTimeDependentDirichletBCs ()
    {
        dolfin::begin (dolfin::DBG, "Setting time dependent Dirichlet's boudary conditions...");
        
        // reset time dependent Dirichlet BCs
        for (auto bcIterator = timeDependentDirichletBCs_.begin (); 
             bcIterator != timeDependentDirichletBCs_.end (); 
             bcIterator++)
        {
            resetTimeDependentDirichletBC (bcIterator);
        }

        dolfin::end ();
    }

    
    
    void MovingTimeDependentProblem::resetTimeDependentDirichletBC 
        (std::map <MovingTimeDependentProblem::TimeDependentDirichletBCKey, 
                   MovingTimeDependentProblem::TimeDependentDirichletBCValue>
                   ::iterator bcIterator)
    {
        // get bc name
        std::string bcName = bcIterator->first;
        
        // get dcp::TimeDependentExpression object and set time equal to t_
        std::shared_ptr<dcp::TimeDependentExpression> condition = std::get<0> (bcIterator->second);
        condition->setTime (t_);
        
        // get dcp::Subdomain object
        std::shared_ptr<dcp::Subdomain> boundary = std::get<1> (bcIterator->second);

        // get component on which to enforce the bc
        int component = std::get<2> (bcIterator->second);
        
        dolfin::begin (dolfin::DBG, 
                       "Resetting time dependent Dirichlet's boundary condition \"%s\"...", 
                       bcName.c_str ());
        
        // remove the bc from timeSteppingProblem_
        removeDirichletBC (bcName);
        
        // remember that component = -1 means that the bc should be enforced on all the function space components
        if (component == -1)
        {
            addDirichletBC (condition, boundary, bcName);
        }
        else
        {
            addDirichletBC (condition, boundary, component, bcName);
        } 
        
        dolfin::end ();
    }

    

    void MovingTimeDependentProblem::setTimeDependentCoefficients ()
    {
        dolfin::begin (dolfin::DBG, "Setting time dependent coefficients...");
        
        for (auto& coefficientPair : timeDependentCoefficients_)
        {
            std::shared_ptr <dcp::TimeDependentExpression> expression (coefficientPair.second);
            std::string coefficientName = std::get<0> (coefficientPair.first);
            std::string coefficientType = std::get<1> (coefficientPair.first);
            dolfin::begin (dolfin::DBG, 
                         "Coefficient: name \"%s\", type \"%s\"", 
                         coefficientName.c_str (),
                         coefficientType.c_str ());
            
            dolfin::log (dolfin::DBG, "Setting time in time dependent expression...");
            expression->setTime (t_);
            
            timeSteppingProblem_->setCoefficient (coefficientType, expression, coefficientName);
            
            dolfin::end ();
        }
        
        dolfin::end ();
    }
    
            
    
    void MovingTimeDependentProblem::printFinishedWarning ()
    {
        dolfin::warning ("No time iteration performed in solve() function. End time already reached.");
    }
    


    void MovingTimeDependentProblem::storeSolution (const dolfin::Function& solution, 
                                              const int& timeStep, 
                                              const int& storeInterval)
    {
        if (storeInterval > 0 && timeStep % storeInterval == 0)
        {
            dolfin::log (dolfin::DBG, "Saving time stepping problem solution in solutions vector...");
            solution_.push_back (std::make_pair (t_, solution));
        }
    } 
    


    void MovingTimeDependentProblem::storeLastStepSolution (const dolfin::Function& solution, 
                                                      const int& timeStep, 
                                                      const int& storeInterval)
    {
        if (!(storeInterval > 0 && timeStep % storeInterval == 0))
        {
            dolfin::log (dolfin::DBG, "Saving last time step solution in solutions vector...");
            solution_.push_back (std::make_pair (t_, solution));
        }
    }

    

    void MovingTimeDependentProblem::plotSolution (dolfin::Function& solution, 
                                             const int& timeStep, 
                                             const int& plotInterval, 
                                             const int& plotComponent,
                                             const bool& pause)
    {
        if (plotInterval > 0 && timeStep % plotInterval == 0)
        {
            dolfin::begin (dolfin::DBG, "Plotting...");

            // auxiliary variable, to enhance readability
            std::shared_ptr<dolfin::Function> functionToPlot;
            
            std::string plotTitle = parameters ["plot_title"];

            // if plotTitle is not empty, we need to prepend ", " so that the plot title is readable
            if (!plotTitle.empty ())
            {
                plotTitle = ", " + plotTitle;
            }

            // get right function to plot
            if (plotComponent == -1)
            {
                functionToPlot = dolfin::reference_to_no_delete_pointer (solution);
                dolfin::log (dolfin::DBG, "Plotting time stepping problem solution, all components...");
            }
            else
            {
                functionToPlot = dolfin::reference_to_no_delete_pointer (solution [plotComponent]);
                dolfin::log (dolfin::DBG, "Plotting time stepping problem solution, component %d...", plotComponent);
            }

            // actual plotting
            if (solutionPlotter_ == nullptr)
            {
                dolfin::log (dolfin::DBG, "Plotting in new dolfin::VTKPlotter object...");
                solutionPlotter_ = dolfin::plot (functionToPlot, "Time = " + std::to_string (t_) + plotTitle);
            }
            else if (! solutionPlotter_ -> is_compatible (functionToPlot))
            {
                dolfin::log (dolfin::DBG, "Existing plotter is not compatible with object to be plotted.");
                dolfin::log (dolfin::DBG, "Creating new dolfin::VTKPlotter object...");
                solutionPlotter_ = dolfin::plot (functionToPlot, "Time = " + std::to_string (t_) + plotTitle);
            }
            else 
            {
                solutionPlotter_ -> parameters ["title"] = std::string ("Time = " + std::to_string (t_) + plotTitle);
                solutionPlotter_ -> plot (functionToPlot);
            }

            if (pause)
            {
                dolfin::interactive ();
            }

            dolfin::end ();
        }
    }
}
