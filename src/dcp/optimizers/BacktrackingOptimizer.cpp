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
#include <dcp/utils/dotproductforms.h>

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

        dolfin::end (); // Creating BacktrackingOptimizer object
        
        dolfin::log (dolfin::DBG, "BacktrackingOptimizer object created");
        
        dolfin::end ();
    }
    

    // ***************************************** //
    // ********** PRIVATE MEMBERS ************** //
    // ***************************************** //
    void BacktrackingOptimizer::print_ (std::ostream& OUTSTREAM, 
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
            OUTSTREAM << "  "
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
            OUTSTREAM << "  "
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
            OUTSTREAM << "  "
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
    


    bool BacktrackingOptimizer::openOutputFile_ (std::ofstream& outfile) const
    {
        std::string outputFileName = this->parameters ["output_file_name"];
        std::string convergenceCriterion = this->parameters ["convergence_criterion"];
        
        if (!outputFileName.empty ())
        {
            outfile.open (outputFileName);
            if (outfile.fail ())
            {
                dolfin::dolfin_error ("dcp: BacktrackingOptimizer.cpp",
                                      "openOutputFile_",
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
                                      "openOutputFile_",
                                      "Unknown convergence criterion \"%s\"", 
                                      convergenceCriterion.c_str ());
                
            }
            return true;
        }
        
        return false;
    }
}
