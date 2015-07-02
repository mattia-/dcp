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

#include <dcp/differential_problems/GenericProblem.h>
#include <dolfin/log/dolfin_log.h>
#include <map>
#include <string>
#include <utility>

namespace dcp
{
    /************************* CONSTRUCTORS ********************/
    GenericProblem::GenericProblem (const std::shared_ptr<dolfin::FunctionSpace> functionSpace) : 
        parameters ("differential_problem_parameters"),
        functionSpace_ (functionSpace),
        dirichletBCs_ (),
        solution_ (),
        stashedSolution_ (*functionSpace_),
        dirichletBCsCounter_ (0),
        solutionPlotter_ ()
    { 
        dolfin::begin (dolfin::DBG, "Building GenericProblem...");
        
        dolfin::log (dolfin::DBG, "Setting up parameters...");
        parameters.add ("plot_component", -1);
        parameters.add ("plot_title", "Solution");
        parameters.add ("clone_method", "shallow_clone");
            
        dolfin::end ();
        
        dolfin::log (dolfin::DBG, "GenericProblem object created");
    }



    GenericProblem::GenericProblem (const dolfin::FunctionSpace& functionSpace) : 
        parameters ("differential_problem_parameters"),
        functionSpace_ (new dolfin::FunctionSpace (functionSpace)),
        dirichletBCs_ (),
        solution_ (),
        stashedSolution_ (functionSpace_),
        dirichletBCsCounter_ (0),
        solutionPlotter_ ()
    { 
        dolfin::begin (dolfin::DBG, "Building GenericProblem...");
        
        dolfin::log (dolfin::DBG, "Setting up parameters...");
        parameters.add ("plot_component", -1);
        parameters.add ("plot_title", "Solution");
        parameters.add ("clone_method", "shallow_clone");
            
        dolfin::end ();
        
        dolfin::log (dolfin::DBG, "GenericProblem object created");
    }



    GenericProblem::GenericProblem (dolfin::FunctionSpace&& functionSpace) : 
        parameters ("differential_problem_parameters"),
        functionSpace_ (new dolfin::FunctionSpace (std::move (functionSpace))),
        dirichletBCs_ (),
        solution_ (),
        stashedSolution_ (functionSpace_),
        dirichletBCsCounter_ (0),
        solutionPlotter_ ()
    { 
        dolfin::begin (dolfin::DBG, "Building GenericProblem...");
        
        dolfin::log (dolfin::DBG, "Setting up parameters...");
        parameters.add ("plot_component", -1);
        parameters.add ("plot_title", "Solution");
        parameters.add ("clone_method", "shallow_clone");
            
        dolfin::end ();
        
        dolfin::log (dolfin::DBG, "GenericProblem object created");
    }



    /********************** GETTERS ***********************/
    std::shared_ptr<const dolfin::Mesh> GenericProblem::mesh () const
    {
        return functionSpace_ -> mesh ();      
    }



    std::shared_ptr<dolfin::FunctionSpace> GenericProblem::functionSpace () const
    {
        return functionSpace_;
    }



    const dolfin::DirichletBC& GenericProblem::dirichletBC (const std::string& bcName) const
    {
        auto bcIterator = dirichletBCs_.find (bcName);
        if (bcIterator == dirichletBCs_.end ())
        {
            dolfin::dolfin_error ("dcp: GenericProblem.cpp",
                                  "dirichletBC",
                                  "Cannot find dirichletBC with name \"%s\" in map", 
                                  bcName.c_str ());
        }
        return bcIterator -> second;
    }



    const std::map<std::string, dolfin::DirichletBC>& GenericProblem::dirichletBCs () const
    {
        return dirichletBCs_;
    }



    const dolfin::Function& GenericProblem::solution (const std::string& solutionType) const
    {
        if (solutionType == "default")
        {
            return solution_.back ().second;
        }
        else if (solutionType == "stashed")
        {
            return stashedSolution_;
        }
        else
        {
            dolfin::dolfin_error ("dcp: GenericProblem.cpp",
                                  "solution",
                                  "Unkown solution type \"%s\" requested", solutionType.c_str ());
            return stashedSolution_; // just to suppress the compilation warning, dolfin_error will cause the program
                                     // to exit anyway
        }
    }



    /********************** SETTERS ***********************/
    bool GenericProblem::addDirichletBC (const dolfin::GenericFunction& condition, 
                                         const dolfin::SubDomain& boundary,
                                         std::string bcName)
    {
        return addDirichletBC (dolfin::reference_to_no_delete_pointer (condition), 
                               dolfin::reference_to_no_delete_pointer (boundary),
                               bcName); 
    }
    


    bool GenericProblem::addDirichletBC (const dolfin::GenericFunction& condition, 
                                         const dolfin::SubDomain& boundary, 
                                         const std::size_t& component,
                                         std::string bcName)
    {
        return addDirichletBC (dolfin::reference_to_no_delete_pointer (condition), 
                               dolfin::reference_to_no_delete_pointer (boundary),
                               component,
                               bcName); 
    }
    


    bool GenericProblem::addDirichletBC (std::shared_ptr<const dolfin::GenericFunction> condition, 
                                         std::shared_ptr<const dolfin::SubDomain> boundary,
                                         std::string bcName)
    {
        if (bcName.empty ())
        {
            bcName = "dirichlet_condition_" + std::to_string (dirichletBCsCounter_);
            dirichletBCsCounter_++;
        }
        
        dolfin::log (dolfin::DBG, 
                     "Adding dirichlet boundary condition to boundary conditions map with name \"%s\"...",
                     bcName.c_str ());
        
        auto result = dirichletBCs_.emplace (bcName, dolfin::DirichletBC (functionSpace_, condition, boundary));
        
        if (result.second == false)
        {
            dolfin::warning ("DirichletBC object not inserted because key \"%s\" already in map",
                             bcName.c_str ());
        }
        
        return result.second;
    }

    

    bool GenericProblem::addDirichletBC (std::shared_ptr<const dolfin::GenericFunction> condition, 
                                         std::shared_ptr<const dolfin::SubDomain> boundary,
                                         const std::size_t& component,
                                         std::string bcName)
    {
        if (bcName.empty ())
        {
            bcName = "dirichlet_condition_" + std::to_string (dirichletBCsCounter_);
            dirichletBCsCounter_++;
        }
        
        dolfin::log (dolfin::DBG, 
                     "Adding dirichlet boundary condition to boundary conditions map with name \"%s\"...",
                     bcName.c_str ());
        
        auto result = dirichletBCs_.emplace 
            (bcName, dolfin::DirichletBC ((*functionSpace_) [component], condition, boundary));
        
        if (result.second == false)
        {
            dolfin::warning ("DirichletBC object not inserted because key \"%s\" already in map",
                             bcName.c_str ());
        }
        
        return result.second;
    }

    

    bool GenericProblem::addDirichletBC (const dolfin::DirichletBC& dirichletCondition, 
                                         std::string bcName)
    {
        if (bcName.empty ())
        {
            bcName = "dirichlet_condition_" + std::to_string (dirichletBCsCounter_);
            dirichletBCsCounter_++;
        }
        
        dolfin::log (dolfin::DBG, 
                     "Adding dirichlet boundary condition to boundary conditions map with name \"%s\"...",
                     bcName.c_str ());
        
        auto result = dirichletBCs_.insert (std::make_pair (bcName, dirichletCondition));
        
        if (result.second == false)
        {
            dolfin::warning ("DirichletBC object not inserted because key \"%s\" already in map",
                             bcName.c_str ());
        }
        
        return result.second;
    }



    bool GenericProblem::addDirichletBC (dolfin::DirichletBC&& dirichletCondition,
                                         std::string bcName)
    {
        if (bcName.empty ())
        {
            bcName = "dirichlet_condition_" + std::to_string (dirichletBCsCounter_);
            dirichletBCsCounter_++;
        }
        
        dolfin::log (dolfin::DBG, 
                     "Adding dirichlet boundary condition to boundary conditions map with name \"%s\"...",
                     bcName.c_str ());
        
        auto result = dirichletBCs_.insert (std::make_pair (bcName, dirichletCondition));
        
        if (result.second == false)
        {
            dolfin::warning ("DirichletBC object not inserted because key \"%s\" already in map", 
                             bcName.c_str ());
        }
        
        return result.second;
    }



    bool GenericProblem::removeDirichletBC (const std::string& bcName)
    {
        dolfin::log (dolfin::DBG, 
                     "Removing dirichlet boundary condition \"%s\" from boundary conditions map...", 
                     bcName.c_str ());
        std::size_t nErasedElements = dirichletBCs_.erase (bcName);
        
        if (nErasedElements == 0)
        {
            dolfin::warning ("Dirichlet boundary condition \"%s\" not found in map", 
                             bcName.c_str ());
        }
        
        return nErasedElements == 1? true : false;
    }
    


    void GenericProblem::update ()
    {

    }
    
    

    /********************** METHODS ***********************/
    void GenericProblem::plotSolution (const std::string& plotType)
    {
        // check if plotType is known
        if (plotType != "all")
        {
            dolfin::warning ("Uknown plot type \"%s\". No plot performed", plotType.c_str ());
            return;
        }
        
        dolfin::begin (dolfin::DBG, "Plotting...");
        int plotComponent = parameters ["plot_component"];
        
        // auxiliary variable, to enhance readability
        std::shared_ptr<dolfin::Function> functionToPlot;
        
        // get right function to plot
        if (plotComponent == -1)
        {
            functionToPlot = dolfin::reference_to_no_delete_pointer (solution_.back ().second);
            dolfin::log (dolfin::DBG, "Plotting problem solution, all components...");
        }
        else
        {
            functionToPlot = dolfin::reference_to_no_delete_pointer (solution_.back ().second [plotComponent]);
            dolfin::log (dolfin::DBG, "Plotting problem solution, component %d...", plotComponent);
        }
        
        // actual plotting
        if (solutionPlotter_ == nullptr)
        {
            dolfin::log (dolfin::DBG, "Plotting in new dolfin::VTKPlotter object...");
            solutionPlotter_ = dolfin::plot (functionToPlot, parameters ["plot_title"]);
        }
        else if (! solutionPlotter_ -> is_compatible (functionToPlot))
        {
            dolfin::log (dolfin::DBG, "Existing plotter is not compatible with object to be plotted.");
            dolfin::log (dolfin::DBG, "Creating new dolfin::VTKPlotter object...");
            solutionPlotter_ = dolfin::plot (functionToPlot, parameters ["plot_title"]);
        }
        else 
        {
            solutionPlotter_ -> parameters ["title"] = parameters ["plot_title"];
            solutionPlotter_ -> plot (functionToPlot);
        }
        
        dolfin::end ();
    }
    


    void GenericProblem::applyStashedSolution ()
    {
        solution_.back ().second = stashedSolution_;
        stashedSolution_ = dolfin::Function (*functionSpace_);
    }
}
