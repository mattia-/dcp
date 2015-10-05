/* 
 *  Copyright (C) 2015, Ivan Fumagalli, ivan.fumagalli.if@gmail.com
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

#define EVITAPLOT
//#define MULTISTEP

//#define COMPUTESTRESS

//#include <dcp/differential_problems/MovingTimeDependentProblem.h>
#include "MovingTimeDependentProblem.h"
#include <dolfin/log/dolfin_log.h>
#include <regex>
//#include "utilities.h"
// per computeFreeSurfaceStress
#ifdef COMPUTESTRESS
#include "computeFreeSurfaceStress_onlyTP.h"
#include "computeFreeSurfaceStress_noTP.h"
#endif

namespace Ivan
{
    /******************* CONSTRUCTORS *******************/
    MovingTimeDependentProblem::MovingTimeDependentProblem 
//        (const std::shared_ptr<geometry::MeshManager<dolfin::ALE,dolfin::FunctionSpace> > meshManager,
        (const std::shared_ptr<MeshManager<> > meshManager,
			 	 const std::shared_ptr<dcp::AbstractProblem> timeSteppingProblem,
         const double& startTime,
         const double& dt,
         const double& endTime,
         std::initializer_list<std::string> dtCoefficientTypes,
         std::initializer_list<std::string> previousSolutionCoefficientTypes,
         std::initializer_list<std::string> wCoefficientTypes,
         const unsigned int& nTimeSchemeSteps)
        : 
            TimeDependentProblem (timeSteppingProblem,
                                  startTime, dt, endTime,
                                  dtCoefficientTypes,
                                  previousSolutionCoefficientTypes,
                                  nTimeSchemeSteps),
            meshManager_ (meshManager) //std::shared_ptr<dolfin::ALE>(new dolfin::ALE()),timeSteppingProblem->functionSpace())
    { 
        dolfin::begin (dolfin::DBG, "Building MovingTimeDependentProblem...");
        
        dolfin::Parameters wCoefficientTypesParameter ("w_coefficient_types");
        for (auto& i : wCoefficientTypes)
        {
            wCoefficientTypesParameter.add<std::string> (i);
        }
        parameters.add (wCoefficientTypesParameter);

//        parameters.add ("time_stepping_solution_component", -1);
        parameters["time_stepping_solution_component"] = 0;
          // TODO rimettere default a -1 e capire come passarlo da fuori
        
        dolfin::end ();
        
        dolfin::log (dolfin::DBG, "MovingTimeDependentProblem object created");
    }



    /******************* GETTERS *******************/
    std::shared_ptr<const dolfin::Mesh> MovingTimeDependentProblem::mesh () const
    {
//        return timeSteppingProblem_ -> mesh ();
        return meshManager_->mesh ();
    }



    std::shared_ptr<dolfin::FunctionSpace> MovingTimeDependentProblem::functionSpace () const
    {
        return timeSteppingProblem_ -> functionSpace ();
    }



    /******************* METHODS *******************/
    
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
        std::vector<std::string> wCoefficientTypes;
        parameters ("w_coefficient_types").get_parameter_keys (wCoefficientTypes);
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
        dolfin::Function oldSolution (tmpSolution);
        dolfin::Function displacement (tmpSolution);
        
        std::shared_ptr<dolfin::VTKPlotter> plotter;
        
// per computeFreeSurfaceStress
#ifdef COMPUTESTRESS
//lasciarlo qui?? l'importante è che il marking avvenga solo una volta (quindi prima del ciclo temporale)
// è per calcolare la tensione superficiale
    dolfin::FacetFunction<std::size_t> meshFacets (meshManager_->mesh());
    meshFacets.set_all (0);
    TopBd freeSurface;
    freeSurface.mark (meshFacets, 1);
    LateralBd wallBoundary;
    wallBoundary.mark (meshFacets, 3);
TriplePointLeft triplePointLeft;
triplePointLeft.mark(meshFacets,4);
TriplePointRight triplePointRight;
triplePointRight.mark(meshFacets,5);
dolfin::plot(meshFacets,"MovingTimeDepPb meshFacets");//dolfin::interactive();
#endif
											
        // start time loop. The loop variable timeStepFlag is true only if the solve type requested is "step" and
        // the first step has already been peformed
        int timeStep = 0;
        bool timeStepFlag = false;
dolfin::File solutionFile(problemData.savepath+"insideMovingTimeDependentProblem.pvd");
solutionFile << std::pair<dolfin::Function*,double>(&tmpSolution,t_);
dolfin::HDF5File solutionHDF5File(MPI_COMM_WORLD,problemData.savepath+"sol.hdf5","w");
solutionHDF5File.write (tmpSolution,std::to_string(timeStep));
// per computeFreeSurfaceStress
#ifdef COMPUTESTRESS
dolfin::File surfaceStressFile_onlyTP(problemData.savepath+"surfaceTension_onlyTP.pvd");
//std::ofstream surfaceStressFile_onlyTP_csv; surfaceStressFile_onlyTP_csv.open(savepath+"surfaceTension_onlyTP.csv",std::ios::out|std::ios::app);
dolfin::File surfaceStressFile_noTP(problemData.savepath+"surfaceTension_noTP.pvd");
//dolfin::File provaFile(savepath+"provaFile.pvd");
#endif
dolfin::File displacementFile(problemData.savepath+"displacement.pvd");
const dolfin::Function * w (meshManager_->displacement().get());
w = meshManager_->computeDisplacement(displacement[0],"normal",dt).get();
displacementFile << std::make_pair(w,t_);

dolfin::Constant zero (0,0);
        while (!isFinished () && timeStepFlag == false)
        {
            timeStep++;
            
            //meshManager_->moveMesh(tmpSolution[0],"normal",dt);
            meshManager_->moveMesh();

            for (auto& i : wCoefficientTypes)
            {
                timeSteppingProblem_->setCoefficient (i, dolfin::reference_to_no_delete_pointer (*w), "w");
//std::cerr << i << std::endl;
            } 
        /*    timeSteppingProblem_->setCoefficient ("residual_form", dolfin::reference_to_no_delete_pointer (w), "w");
            timeSteppingProblem_->setCoefficient ("jacobian_form", dolfin::reference_to_no_delete_pointer (w), "w");*/
            
            if (oneStepRequested == false)
            {
                dolfin::begin (dolfin::INFO, "===== Timestep %d =====", timeStep);
            }
            
            advanceTime ();
            
            setTimeDependentCoefficients ();
            
            setTimeDependentDirichletBCs ();
            
            dolfin::begin (dolfin::INFO, "Solving time stepping problem...");
// per computeFreeSurfaceStress
#ifdef COMPUTESTRESS
//provaFile << timeSteppingProblem_->solution ();
#endif
            timeSteppingProblem_->solve ();
            
            oldSolution = tmpSolution;
            tmpSolution = timeSteppingProblem_->solution ();
solutionFile << std::pair<dolfin::Function*,double>(&tmpSolution,t_);
solutionHDF5File.write (tmpSolution,std::to_string(timeStep));

// per computeFreeSurfaceStress
#ifdef COMPUTESTRESS
// Assembling and storing the surface tension term
computeFreeSurfaceStress_noTP::FunctionSpace stressSpace_noTP(meshManager_->mesh());
computeFreeSurfaceStress_onlyTP::FunctionSpace stressSpace(meshManager_->mesh());
dolfin::Constant stressGamma (problemData.gamma);
dolfin::Constant stressDt (problemData.dt);
dolfin::Constant beta (problemData.beta);
dolfin::Constant cosThetaS (cos ( problemData.thetaS * 3.14159265/180.0 ));
dolfin::Constant t_partialOmega(0,1);
// CI FOSSE residualForm() COME METODO PER ESTRARRE I COEFFICIENTI...  stressForm.set_coefficient ("gamma",timeSteppingProblem_->residualForm().coefficient("gamma"));
computeFreeSurfaceStress_noTP::LinearForm stressForm_noTP(stressSpace_noTP,tmpSolution,stressDt,stressGamma,beta,*(new WallVelocity));
computeFreeSurfaceStress_onlyTP::LinearForm stressForm_onlyTP(stressSpace,stressDt,stressGamma,cosThetaS,t_partialOmega);
stressForm_noTP.set_exterior_facet_domains (dolfin::reference_to_no_delete_pointer(meshFacets));
stressForm_onlyTP.set_exterior_facet_domains (dolfin::reference_to_no_delete_pointer(meshFacets));
//stressForm.set_vertex_domains (dolfin::reference_to_no_delete_pointer(meshPoints));
dolfin::Vector b;
dolfin::Function bFun(stressSpace);
assemble(b, stressForm_onlyTP);            
* bFun.vector() = b;
surfaceStressFile_onlyTP << std::pair<dolfin::Function*,double>(&bFun,t_);
//for (std::size_t i=0; i!=b.size(); ++i) surfaceStressFile_onlyTP_csv << b[i] << ","; surfaceStressFile_onlyTP_csv << std::endl;
dolfin::Function bFun_noTP(stressSpace_noTP);
assemble(* bFun_noTP.vector(), stressForm_noTP);
surfaceStressFile_noTP << std::pair<dolfin::Function*,double>(&bFun_noTP,t_);
std::cerr << std::endl;
#endif

            // save solution in solution_ according to time step and store interval.
            storeSolution (tmpSolution, timeStep, storeInterval);
            
#ifndef EVITAPLOT
            // plot solution according to time step and plot interval
 //           plotSolution (tmpSolution, timeStep, plotInterval, plotComponent, pause);
            plotSolution (tmpSolution, timeStep, 1, 0, false);//dolfin::interactive();
//            plotSolution (tmpSolution, timeStep, 1, 1, false);//dolfin::interactive();
            dolfin::plot (*meshManager_->mesh(),"inside mesh");//dolfin::interactive();
#endif

            if (oneStepRequested == true)
            {
                timeStepFlag = true;
            }
             
            // dolfin::end matching dolfin::begin inside if statement
            if (oneStepRequested == false)
            {
                dolfin::end ();
            }
        
            //w =
            //meshManager_->computeDisplacement(tmpSolution[0],"normal",dt).get();
//for (unsigned i(0); i<5; ++i) std::cerr << (*tmpSolution.vector())[i] << ' '; std::cerr << std::endl;
//for (unsigned i(0); i<5; ++i) std::cerr << (*oldSolution.vector())[i] << ' '; std::cerr << std::endl;
            displacement = tmpSolution;
#ifdef MULTISTEP
            *displacement.vector() *= 2;
            *displacement.vector() -= * oldSolution.vector();
#endif
//for (unsigned i(0); i<5; ++i) std::cerr << (*displacement.vector())[i] << ' '; std::cerr << std::endl;
//for (unsigned i(0); i<5; ++i) std::cerr << (*tmpSolution.vector())[i] << ' '; std::cerr << std::endl;
            w = meshManager_->computeDisplacement(displacement[0],"normal",dt).get();
            displacementFile << std::make_pair(w,t_);
#ifndef EVITAPLOT
            dolfin::plot (*w, "w"); //dolfin::interactive();
#endif

            dolfin::end ();
        }
        
        // At this point, we just need to make sure that the solution on the last iteration was saved even though
        // timeStep % storeInterval != 0 (but it must not be saved twice!)
        storeLastStepSolution (tmpSolution, timeStep, storeInterval);

// per computeFreeSurfaceStress
#ifdef COMPUTESTRESS
//surfaceStressFile_onlyTP_csv.close();
#endif
solutionHDF5File.close();
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
                new Ivan::MovingTimeDependentProblem (this->meshManager_, this->timeSteppingProblem_, startTime_, dt_, endTime_, {}, {}, {});
            clonedProblem->timeDependentDirichletBCs_ = this->timeDependentDirichletBCs_;
        }
        else if (cloneMethod == "deep_clone")
        {
            // note that we pass an empty initializer_list to the constructor as dtCoefficientTypes and 
            // previousSolutionCoefficientTypes, because they will be copied when the parameters are copied anyway
            clonedProblem = 
                new Ivan::MovingTimeDependentProblem (this->meshManager_, this->timeSteppingProblem_, startTime_, dt_, endTime_, {}, {}, {});
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
    
    void MovingTimeDependentProblem::initializeMesh (dolfin::Expression & displacement)
    {
        dolfin::Function displ (solution_.back().second.function_space());
        displ = displacement;
        meshManager_->moveMesh (displ[0]);
    }


} //end of namespace
