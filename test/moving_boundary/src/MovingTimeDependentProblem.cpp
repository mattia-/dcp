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
//#define PARAB
//#define TUTTELECOMP
#define POSTPROCESSOR

//#include <dcp/differential_problems/MovingTimeDependentProblem.h>
#include "MovingTimeDependentProblem.h"
#include <dolfin/log/dolfin_log.h>
#include <regex>
//#include "utilities.h"

// TODO toglierle: servono solo per il preassemble
//#include "heateq.h"
//#include "heateqOld.h"
#include "myNavierstokesTimeCurvLinear.h"
#include "myNavierstokesTimeCurvLinearPreviousDomain.h"
#include "computeFreeSurfaceStress_onlyTP.h"

#undef EVITAPLOT

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
            meshManager_ (meshManager), //std::shared_ptr<dolfin::ALE>(new dolfin::ALE()),timeSteppingProblem->functionSpace())
            postProcessor_ (* this)
    { 
        dolfin::begin (dolfin::DBG, "Building MovingTimeDependentProblem...");
        
        dolfin::log (dolfin::DBG, "Setting up parameters...");
        parameters.add ("mesh_displacement_name", "w");
        dolfin::Parameters wCoefficientTypesParameter ("mesh_displacement_coefficient_types");
        for (auto& i : wCoefficientTypes)
        {
            wCoefficientTypesParameter.add<std::string> (i);
        }
        parameters.add (wCoefficientTypesParameter);

#ifdef TUTTELECOMP
        parameters["time_stepping_solution_component"] = -1;
        parameters["plot_component"] = -1;
#else
        parameters["time_stepping_solution_component"] = 0;
#endif
          // TODO rimettere default a -1 e capire come passarlo da fuori
        
        dolfin::end ();
        
        dolfin::log (dolfin::DBG, "MovingTimeDependentProblem object created");

#ifdef MULTISTEP
std::cerr << " ----- extrapolation of ALE velocity activated -----" << std::endl;
#endif
    }


    MovingTimeDependentProblem::~MovingTimeDependentProblem ()
    {
#ifdef MULTISTEP
std::cerr << " ----- extrapolation of ALE velocity was used -----" << std::endl;
#endif
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



    MeshManager<> & MovingTimeDependentProblem::meshManager () const
    {
        return * meshManager_;
    }


    /******************* SETTERS *******************/
    void MovingTimeDependentProblem::initializeMesh (dolfin::Expression & displacement)
    {
        dolfin::Function displ (solution_.back().second.function_space());
        displ = displacement;
        meshManager_->moveMesh (displ[0], "init", 1);
    }


    /******************* METHODS *******************/
    
    void MovingTimeDependentProblem::solve (const std::string& type) 
    {

        std::string solFileName;

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
        parameters ("mesh_displacement_coefficient_types").get_parameter_keys (wCoefficientTypes);
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
        
        
        // functions used to step through the time loop
        dolfin::Function tmpSolution = solution_.back ().second;
        dolfin::Function oldSolution (tmpSolution);
        dolfin::Function displacement (tmpSolution);
		const dolfin::Function * w (meshManager_->displacement().get());  // non-const pointer to const dolfin::Function
		//w = meshManager_->computeDisplacement(displacement,"normal",dt).get();
		w = meshManager_->computeDisplacement(displacement[0],"normal",dt).get();
      //TODO questa componente va passata da fuori
        
        std::shared_ptr<dolfin::VTKPlotter> plotter;

        // start time loop. The loop variable timeStepFlag is true only if the solve type requested is "step" and
        // the first step has already been peformed
        int timeStep = 0;
        bool timeStepFlag = false;
#ifdef POSTPROCESSOR
postProcessor_ (timeStep);
#endif

if (timeStep % problemData.fileFrequency <= 2)
{
    if (problemData.saveDataInFile)
  		saveDataInFile (tmpSolution, *w, timeStep);
    if (problemData.savePvd)
    {
      solFileName = problemData.savepath+"solution"+std::to_string(timeStep)+".pvd";
      dolfin::File solFile (solFileName);
      solFile << tmpSolution;
    }
}

        while (!isFinished () && timeStepFlag == false)
			//TODO non mi piace questa flag: magari gia' e' stata tolta in un'altra branch?
        {

            timeStep++;
            
            if (oneStepRequested == false)
            {
                dolfin::begin (dolfin::INFO, "===== Timestep %d =====", timeStep);
            }
            
            advanceTime ();

if (timeStep % problemData.fileFrequency <= 2)
{
#ifdef POSTPROCESSOR
postProcessor_.onOldDomain (timeStep);
#endif
}

            //std::static_pointer_cast<
            //  Ivan::MovingLinearProblem <heateq::FunctionSpace, heateq::FunctionSpace, heateq::BilinearForm, heateq::LinearForm, heateqOld::LinearForm, heateqOld::LinearForm> >
           // std::static_pointer_cast<Ivan::MovingLinearProblem < myNavierstokesTimeCurvLinear::FunctionSpace, computeFreeSurfaceStress_onlyTP::FunctionSpace,
           //           myNavierstokesTimeCurvLinear::BilinearForm, myNavierstokesTimeCurvLinear::LinearForm,
           //           myNavierstokesTimeCurvLinearPreviousDomain::LinearForm, computeFreeSurfaceStress_onlyTP::LinearForm > >
           std::dynamic_pointer_cast<Ivan::MovingAbstractProblem>
                (timeSteppingProblem_)->preassemble ();
			// TODO considerare la presenza di un preassemble in AbstractProblem, per evitare casting
			//      o magari usare la solve() passando "preassemble" come stringa
            //meshManager_->moveMesh(tmpSolution[0],"normal",dt);
            meshManager_->moveMesh();
/*std::cerr << "displ : ";
for (std::size_t i (0); i!=w->vector()->size(); ++i)
  std::cerr << std::setprecision(10) << (* w->vector())[i] << ", ";
std::cerr << std::endl;*/
            //this->adapt ();
            
            setTimeDependentCoefficients ();
            
            setTimeDependentDirichletBCs ();
            
//            std::static_pointer_cast<
//              Ivan::MovingLinearProblem <heateq::FunctionSpace, heateq::FunctionSpace, heateq::BilinearForm, heateq::LinearForm, heateqOld::LinearForm, heateqOld::LinearForm> >
//                (timeSteppingProblem_)->preassemble ();
//			// TODO considerare la presenza di un preassemble in AbstractProblem, per evitare casting
//            //meshManager_->moveMesh(tmpSolution[0],"normal",dt);
//            meshManager_->moveMesh();

            for (auto& i : wCoefficientTypes)
            {
                timeSteppingProblem_->setCoefficient (i, dolfin::reference_to_no_delete_pointer (*w), "w");
            } 
            
            dolfin::begin (dolfin::INFO, "Solving time stepping problem...");

if (timeStep % problemData.fileFrequency <= 2)
            timeSteppingProblem_->solve ("save_residual");
else
            timeSteppingProblem_->solve ();
            
			dolfin::end ();

            oldSolution = tmpSolution;
            tmpSolution = timeSteppingProblem_->solution ();
            // save solution in solution_ according to time step and store interval.
            storeSolution (tmpSolution, timeStep, storeInterval);
            
#ifndef EVITAPLOT
            // plot solution according to time step and plot interval
 //           plotSolution (tmpSolution, timeStep, plotInterval, plotComponent, pause);
  #ifdef TUTTELECOMP
            plotSolution (tmpSolution, timeStep, 1, -1, false);//dolfin::interactive();
  #else
            plotSolution (tmpSolution, timeStep, 1, 0, false);//dolfin::interactive();
  #endif
//            plotSolution (tmpSolution, timeStep, 1, 1, false);//dolfin::interactive();
            dolfin::plot (*meshManager_->mesh(),"inside mesh");//dolfin::interactive();
#endif

if (timeStep % problemData.fileFrequency <= 2)
{
#ifdef POSTPROCESSOR
postProcessor_ (timeStep);
#endif
}

            displacement = tmpSolution;
#ifdef MULTISTEP
            *displacement.vector() *= 2;
            *displacement.vector() -= * oldSolution.vector();
#endif
/*#ifdef PARAB
            dolfin::Function displDir (* postProcessor_.uDirSize3_.function_space());
#ifdef TUTTELECOMP
            dolfin::assign (dolfin::reference_to_no_delete_pointer (displDir[0]), dolfin::reference_to_no_delete_pointer (displacement));
#else
            dolfin::assign (dolfin::reference_to_no_delete_pointer (displDir[0]), dolfin::reference_to_no_delete_pointer (displacement[0]));
#endif
            * displDir.vector() += * postProcessor_.uDirSize3_.vector();
            w = meshManager_->computeDisplacement(displDir[0],"normal",dt).get();
#else
            w = meshManager_->computeDisplacement(displacement[0],"normal",dt).get();
#endif
*/
            //w = meshManager_->computeDisplacement(displacement,"normal",dt).get();
            w = meshManager_->computeDisplacement(displacement[0],"normal",dt).get();
              //TODO questa componente va passata da fuori
#ifndef EVITAPLOT
            dolfin::plot (*w, "w"); //dolfin::interactive();
#endif

if (timeStep % problemData.fileFrequency <= 2)
{
    if (problemData.saveDataInFile)
  		saveDataInFile (tmpSolution, *w, timeStep);
    if (problemData.savePvd)
    {
      solFileName = problemData.savepath+"solution"+std::to_string(timeStep)+".pvd";
      dolfin::File solFile (solFileName);
      solFile << tmpSolution;
    }
}

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
	
	
	/******************* PROTECTED METHODS *******************/
    void MovingTimeDependentProblem::saveDataInFile (const dolfin::Function& tmpSolution,
													 const dolfin::Function& w,
													 const int& timeStep) const
	{
print2csv (tmpSolution, problemData.savepath+"sol_p.csv."+std::to_string(timeStep),
           meshManager_->orderedPressDofs().begin(),
           meshManager_->orderedPressDofs().end(),
           tmpSolution.function_space()->dofmap()->tabulate_all_coordinates (* meshManager_->mesh())
           );
print2csv (tmpSolution, problemData.savepath+"sol_vel.csv."+std::to_string(timeStep),
           meshManager_->orderedUXDofs().begin(),
           meshManager_->orderedUYDofs().begin(),
           meshManager_->orderedUXDofs().end(),
           tmpSolution.function_space()->dofmap()->tabulate_all_coordinates (* meshManager_->mesh())
           );
print2csv (w, problemData.savepath+"displacement.csv."+std::to_string(timeStep),
           meshManager_->orderedWXDofs().begin(),
           meshManager_->orderedWYDofs().begin(),
           meshManager_->orderedWXDofs().end(),
           w.function_space()->dofmap()->tabulate_all_coordinates (* meshManager_->mesh())
           );
/*dolfin::HDF5File meshHDF5File(MPI_COMM_WORLD,problemData.savepath+"mesh"+std::to_string(timeStep)+".hdf5","w");
meshHDF5File.write (* tmpSolution.function_space()->mesh(), std::to_string(timeStep));
meshHDF5File.close ();
dolfin::HDF5File solutionHDF5File(MPI_COMM_WORLD,problemData.savepath+"sol"+std::to_string(timeStep)+".hdf5","w");
solutionHDF5File.write (tmpSolution, std::to_string(timeStep));
solutionHDF5File.close ();
dolfin::HDF5File displacementHDF5File(MPI_COMM_WORLD,problemData.savepath+"displ"+std::to_string(timeStep)+".hdf5","w");
displacementHDF5File.write (w, std::to_string(timeStep));
displacementHDF5File.close ();
*/
	}

    void MovingTimeDependentProblem::adapt ()
  {
    //std::static_pointer_cast<
    //  Ivan::MovingLinearProblem <heateq::FunctionSpace, heateq::FunctionSpace, heateq::BilinearForm, heateq::LinearForm, heateqOld::LinearForm, heateqOld::LinearForm> >
          //  std::static_pointer_cast<Ivan::MovingLinearProblem < myNavierstokesTimeCurvLinear::FunctionSpace, computeFreeSurfaceStress_onlyTP::FunctionSpace,
          //            myNavierstokesTimeCurvLinear::BilinearForm, myNavierstokesTimeCurvLinear::LinearForm,
          //            myNavierstokesTimeCurvLinearPreviousDomain::LinearForm, computeFreeSurfaceStress_onlyTP::LinearForm > >
          std::dynamic_pointer_cast<Ivan::MovingAbstractProblem>
            (timeSteppingProblem_)->adapt ();
  }
	
} //end of namespace
