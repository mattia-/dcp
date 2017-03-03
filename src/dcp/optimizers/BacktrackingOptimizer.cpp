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

#include <dcp/optimizers/BacktrackingOptimizer.h>
#include <dcp/utils/dotproductforms.h>
#include <dolfin/parameter/Parameters.h>
#include <dolfin/log/dolfin_log.h>
#include <dolfin/function/Function.h>
#include <dolfin/common/NoDeleter.h>
#include <dolfin/fem/assemble.h>
#include <dolfin/fem/Form.h>
#include <string>
#include <functional>
#include <cmath>
#include <iomanip>

namespace dcp
{
    BacktrackingOptimizer::BacktrackingOptimizer ():
        GenericDescentMethod ()
    {
        dolfin::begin (dolfin::DBG, "Creating BacktrackingOptimizer object...");

        dolfin::log (dolfin::DBG, "Setting up parameters...");
        parameters.add ("descent_method", "backtracking_gradient_method");
        parameters.add ("gradient_norm_tolerance", 1e-6);
        parameters.add ("relative_increment_tolerance", 1e-6);
        parameters.add ("convergence_criterion", "both");
        parameters.add ("c_1", 1e-3);
        parameters.add ("alpha_0", 1.0);
        parameters.add ("rho", 0.5);
        parameters.add ("max_minimization_iterations", 100);
        parameters.add ("max_backtracking_iterations", 20);
        parameters.add ("output_file_name", "");
        parameters.add ("plot_search_direction", false);
        parameters.add ("plot_control_variable", false);
        parameters.add ("pause_internal_plots", false);
        parameters.add ("write_search_direction", false);
        parameters.add ("write_control_variable", false);
        parameters.add ("search_direction_file_name", "search_direction.pvd");
        parameters.add ("control_variable_file_name", "control_variable.pvd");

        dolfin::end (); // Creating BacktrackingOptimizer object

        dolfin::log (dolfin::DBG, "BacktrackingOptimizer object created");
    }


    // ***************************************** //
    // ********** PRIVATE MEMBERS ************** //
    // ***************************************** //
    void BacktrackingOptimizer::printResults_ (std::ostream& outstream,
                                        const int& iteration,
                                        const double& functionalValue,
                                        const double& alpha,
                                        const int& backtrackingIterations,
                                        const double& gradientNorm,
                                        const double& relativeIncrement) const
    {
        dolfin::log (dolfin::DBG, "Printing results to file...");
        std::string convergenceCriterion = this->parameters ["convergence_criterion"];

        if (convergenceCriterion == "both")
        {
            outstream << "  "
                << std::left
                << std::setw (14)
                << iteration
                << std::setw (21)
                << functionalValue
                << std::setw (15)
                << alpha
                << std::setw (28)
                << backtrackingIterations
                << std::setw (18)
                << gradientNorm
                << std::setw (23)
                << relativeIncrement
                << std::endl;
        }
        else if (convergenceCriterion == "increment")
        {
            outstream << "  "
                << std::left
                << std::setw (14)
                << iteration
                << std::setw (21)
                << functionalValue
                << std::setw (15)
                << alpha
                << std::setw (28)
                << backtrackingIterations
                << std::setw (23)
                << relativeIncrement
                << std::endl;
        }
        else if (convergenceCriterion == "gradient")
        {
            outstream << "  "
                << std::left
                << std::setw (14)
                << iteration
                << std::setw (21)
                << functionalValue
                << std::setw (15)
                << alpha
                << std::setw (28)
                << backtrackingIterations
                << std::setw (18)
                << gradientNorm
                << std::endl;
        }
        else
        {
            dolfin::dolfin_error ("dcp: BacktrackingOptimizer.cpp",
                                  "print",
                                  "Unknown convergence criterion \"%s\"",
                                  convergenceCriterion.c_str ());
        }
    }



    bool BacktrackingOptimizer::openResultsFile_ (std::ofstream& outfile) const
    {
        std::string outputFileName = this->parameters ["output_file_name"];
        std::string convergenceCriterion = this->parameters ["convergence_criterion"];

        if (!outputFileName.empty ())
        {
            outfile.open (outputFileName);
            if (outfile.fail ())
            {
                dolfin::dolfin_error ("dcp: BacktrackingOptimizer.cpp",
                                      "openResultsFile_",
                                      "Cannot open output file \"%s\"",
                                      outputFileName.c_str ());
            }

            if (convergenceCriterion == "both")
            {
                outfile << "# Iteration"
                        << "     "
                        << "Functional_value"
                        << "     "
                        << "Alpha"
                        << "          "
                        << "Backtracking_iterations"
                        << "     "
                        << "Gradient_norm"
                        << "     "
                        << "Relative_increment"
                        << std::endl;
            }
            else if (convergenceCriterion == "increment")
            {
                outfile << "# Iteration"
                        << "     "
                        << "Functional_value"
                        << "     "
                        << "Alpha"
                        << "          "
                        << "Backtracking_iterations"
                        << "     "
                        << "Relative_increment"
                        << std::endl;
            }
            else if (convergenceCriterion == "gradient")
            {
                outfile << "# Iteration"
                        << "     "
                        << "Functional_value"
                        << "Alpha"
                        << "          "
                        << "Backtracking_iterations"
                        << "     "
                        << "Gradient_norm"
                        << std::endl;
            }
            else
            {
                dolfin::dolfin_error ("dcp: BacktrackingOptimizer.cpp",
                                      "openResultsFile_",
                                      "Unknown convergence criterion \"%s\"",
                                      convergenceCriterion.c_str ());

            }
            return true;
        }

        return false;
    }



    void BacktrackingOptimizer::setFilenames_
        (const std::vector<const std::shared_ptr<dcp::GenericEquationSystem> > systems,
         std::vector<std::string>& originalFilenames,
         const std::string& action,
         const int& minimizationIteration,
         const int& backtrackingIteration) const
    {
        // clear vector, just for safety
        if (action == "fill_originals")
        {
            originalFilenames.clear ();
        }

        dolfin::begin (dolfin::DBG, "Setting file names for in-loop output to file...");
        // loop through problems in systems and perform the right action according to action
        std::size_t counter = 0;
        for (std::size_t i = 0; i < systems.size (); ++i)
        {
            for (std::size_t j = 0; j < systems[i]->size (); ++j)
            {
                dcp::GenericProblem& problem = (*(systems[i]))[j];

                // regex matching the extension, aka all the characters after the last dot (included)
                std::regex extensionRegex ("(\\.[^.]*$)");

                if (action == "fill_originals")
                {
                    originalFilenames.push_back (problem.parameters["solution_file_name"]);
                }
                else if (action == "set_minimization")
                {
                    // add the minimization iteration number before the extension ($n is the n-th backreference of
                    // the match) in the original filename; remember that names in originalFilenames are stored in
                    // order, so we need to get the counter-th element to have the correct number
                    problem.parameters["solution_file_name"] = std::regex_replace
                        (originalFilenames[counter],
                         extensionRegex,
                         "_minimization_iteration_" + std::to_string (minimizationIteration) + "$1");
                }
                else if (action == "set_backtracking")
                {
                    // add the minimization iteration number and the backtracking iteration number before the
                    // extension ($n is the n-th backreference of the match) in the original filename; remember that
                    // names in originalFilenames are stored in order, so we need to get the counter-th element to
                    // have the correct number
                    problem.parameters["solution_file_name"] = std::regex_replace
                        (originalFilenames[counter],
                         extensionRegex,
                         "_minimization_iteration_" + std::to_string (minimizationIteration)
                         + "_backtracking_iteration_" + std::to_string (backtrackingIteration) + "$1");
                }
                else if (action == "restore_originals")
                {
                    problem.parameters["solution_file_name"] = originalFilenames[counter];
                }
                else
                {
                    dolfin::dolfin_error ("dcp: BacktrackingOptimizer.cpp",
                                          "set filenames",
                                          "Unknown action \"%s\"",
                                          action.c_str ());
                }
            }
            counter++;
        }

        dolfin::end (); // Setting file names for in-loop output to file
    }
}
