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

#include <dcp/problems/GenericProblem.h>
#include <dolfin/log/dolfin_log.h>
#include <map>
#include <string>
#include <utility>
#include <regex>

namespace dcp
{
    /************************* CONSTRUCTORS ********************/
    GenericProblem::GenericProblem (const std::shared_ptr<const dolfin::FunctionSpace> functionSpace) : 
        parameters ("differential_problem_parameters"),
        functionSpace_ (functionSpace),
        dirichletBCs_ (),
        solution_ (),
        stashedSolution_ (functionSpace_),
        dirichletBCsCounter_ (0),
        solutionPlotters_ (),
        solutionFileName_ ("solution.pvd"),
        writeComponents_ (),
        solutionWriters_ ()
    { 
        dolfin::begin (dolfin::DBG, "Building GenericProblem...");
        
        dolfin::log (dolfin::DBG, "Setting up parameters...");
        parameters.add ("solution_file_name", "solution.pvd");
        parameters.add ("plot_components", "-1");
        parameters.add ("plot_title", "Solution");
        parameters.add ("write_components", "-1");
        parameters.add ("clone_method", "shallow_clone");
            
        dolfin::end ();
        
        dolfin::log (dolfin::DBG, "GenericProblem object created");
    }


    /********************** GETTERS ***********************/
    std::shared_ptr<const dolfin::Mesh> GenericProblem::mesh () const
    {
        return functionSpace_ -> mesh ();      
    }



    std::shared_ptr<const dolfin::FunctionSpace> GenericProblem::functionSpace () const
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
    void GenericProblem::applyStashedSolution ()
    {
        solution_.back ().second = stashedSolution_;
        stashedSolution_ = dolfin::Function (functionSpace_);
    }



    void GenericProblem::plotSolution (const std::string& plotType)
    {
        // check if plotType is known
        if (plotType != "default" && plotType != "stashed")
        {
            dolfin::warning ("Unknown plot type \"%s\". No plot performed", plotType.c_str ());
            return;
        }
        
        dolfin::begin (dolfin::DBG, "Plotting...");

        // get vector of plot components
        std::vector<int> plotComponents;
        std::stringstream plotComponentsStream ((std::string (parameters ["plot_components"])));

        // auxiliary variable to push the stream values into the vector
        int component; 
        while (plotComponentsStream >> component)
        {
            plotComponents.push_back (component);
        }
        
        // check if solutionPlotters_ has right size. If not, clear it and reinitialize it with null pointers
        if (solutionPlotters_.size() != plotComponents.size())
        {
            solutionPlotters_.clear();
            solutionPlotters_.resize (plotComponents.size(), nullptr);
        }

        // auxiliary variable, to enhance readability
        std::shared_ptr<dolfin::Function> functionToPlot;
        
        for (std::size_t i = 0; i < plotComponents.size (); ++i)
        {
            int component = plotComponents[i];
 
            // get right function to plot
            if (component == -1)
            {
                if (plotType == "default")
                {
                    functionToPlot = dolfin::reference_to_no_delete_pointer (solution_.back ().second);
                    dolfin::log (dolfin::DBG, "Plotting problem solution, all components...");
                }
                else // aka plotType == "stashed", otherwise we would have exited on the first check
                {
                    functionToPlot = dolfin::reference_to_no_delete_pointer (stashedSolution_);
                    dolfin::log (dolfin::DBG, "Plotting stashed solution, all components...");
                }
            }
            else
            {
                if (plotType == "default")
                {
                    functionToPlot = dolfin::reference_to_no_delete_pointer (solution_.back ().second [component]);
                    dolfin::log (dolfin::DBG, "Plotting problem solution, component %d...", component);
                }
                else // aka plotType == "stashed", otherwise we would have exited on the first check
                {
                    functionToPlot = dolfin::reference_to_no_delete_pointer (stashedSolution_ [component]);
                    dolfin::log (dolfin::DBG, "Plotting stashed solution, component %d...", component);
                }
            }

            std::string plotTitle = parameters ["plot_title"];
            if (component != -1)
            {
                plotTitle += ", component " + std::to_string (component);
            }

            if (plotType == "stashed")
            {
                plotTitle += " (stashed)";
            }

            // actual plotting
            plot_ (solutionPlotters_[i], functionToPlot, plotTitle);
        }
        
        dolfin::end ();
    }
    


    void GenericProblem::writeSolutionToFile (const std::string& writeType)
    {
        // check if writeType is known
        if (writeType != "default" && writeType != "stashed")
        {
            dolfin::warning ("Unknown write type \"%s\". No write performed", writeType.c_str ());
            return;
        }
        
        dolfin::begin (dolfin::DBG, "Saving solution to file...");

        // get vector of write components
        std::vector<int> writeComponents;
        std::stringstream writeComponentsStream ((std::string (parameters ["write_components"])));

        // auxiliary variable to push the stream values into the vector
        int component; 
        while (writeComponentsStream >> component)
        {
            writeComponents.push_back (component);
        }
        
        // check if solutionFileName_ and parameters["solution_file_name"] coincide, if solutionWriters_ has right
        // size and if writeComponents_ contains the same values as writeComponents
        if (solutionFileName_ != std::string (parameters["solution_file_name"])
            ||
            solutionWriters_.size() != writeComponents.size()
            ||
            writeComponents_ != writeComponents)
        {
            solutionFileName_ = std::string (parameters["solution_file_name"]);
            writeComponents_ = writeComponents;

            solutionWriters_.clear ();
            solutionWriters_.resize (writeComponents.size(), nullptr);
        }
        

        // auxiliary variable
        std::shared_ptr<dolfin::Function> functionToWrite;

        for (std::size_t i = 0; i < writeComponents.size (); ++i)
        {
            int component = writeComponents[i];

            if (component == -1)
            {
                if (writeType == "default")
                {
                    functionToWrite = dolfin::reference_to_no_delete_pointer (solution_.back().second);
                }
                else // aka writeType == "stashed", otherwise we would have exited on the first check
                {
                    functionToWrite = dolfin::reference_to_no_delete_pointer (stashedSolution_);
                }
            }

            else
            {
                if (writeType == "default")
                {
                    functionToWrite = dolfin::reference_to_no_delete_pointer (solution_.back().second [component]);
                }
                else // aka writeType == "stashed", otherwise we would have exited on the first check
                {
                    functionToWrite = dolfin::reference_to_no_delete_pointer (stashedSolution_ [component]);
                }
            }

            // file name that keeps track also of the component to be written to file
            std::string filenameWithComponent = solutionFileName_;
            if (component != -1)
            {
                // regex matching the extension, aka all the characters after the last dot (included)
                std::regex extensionRegex ("(\\.[^.]*$)");

                // add the component number before the extension ($n is the n-th backreference of the match)
                filenameWithComponent = std::regex_replace (filenameWithComponent, 
                                                            extensionRegex, 
                                                            "_component" + std::to_string (component) + "$1");
            }

            // actual writing to file
            write_ (solutionWriters_[i], functionToWrite, filenameWithComponent);
        }

        dolfin::end (); // Saving solution to file
    }



    // ---------------------------------------------------------------------------------------------//



    void GenericProblem::plot_ (std::shared_ptr<dolfin::VTKPlotter>& plotter, 
                                const std::shared_ptr<const dolfin::Function> function,
                                const std::string& title)
    {
        if (plotter == nullptr)
        {
            dolfin::log (dolfin::DBG, "Plotting in new dolfin::VTKPlotter object...");
            plotter = dolfin::plot (function, title);
        }
        else if (plotter -> is_compatible (function) == false)
        {
            dolfin::log (dolfin::DBG, "Existing plotter is not compatible with object to be plotted.");
            dolfin::log (dolfin::DBG, "Creating new dolfin::VTKPlotter object...");
            plotter = dolfin::plot (function, title);
        }
        else 
        {
            plotter -> parameters ["title"] = title;
            plotter -> plot (function);
        }
    }



    void GenericProblem::write_ (std::shared_ptr<dolfin::File>& writer, 
                                 const std::shared_ptr<const dolfin::Function> function,
                                 const std::string& filename)
    {
        if (writer == nullptr)
        {
            writer.reset (new dolfin::File (filename));
        }

        (*writer) << (*function);
    }
}
