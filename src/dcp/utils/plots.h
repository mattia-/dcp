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

#ifndef SRC_UTILS_PLOTS_H_INCLUDE_GUARD
#define SRC_UTILS_PLOTS_H_INCLUDE_GUARD

#include <dcp/functions/TimeDependentFunction.h>

// define function with the same name for plots in the dolfin namespace for compatibility
namespace dolfin
{
    //! Dolfin plotter [1]
    /*!
     *  \param function the function to be plotted
     *  \param title the title
     *  \param mode sets the plot mode. Possible values are \c "pause", which pauses the execution after the plot at
     *  each timestep, and \c "continue", which does not. The choice of a string over a simple boolean switch is so that
     *  the signature of the method is exactly the same as for normal dolfin plot functions
     */
    void plot (dcp::TimeDependentFunction& function, std::string title = "", std::string mode = "pause");

    //! Dolfin plotter [2]
    /*!
     *  \param function the function to be plotted
     *  \param title the title
     *  \param mode sets the plot mode. Possible values are \c "pause", which pauses the execution after the plot at
     *  each timestep, and \c "continue", which does not. The choice of a string over a simple boolean switch is so that
     *  the signature of the method is exactly the same as for normal dolfin plot functions
     */
    void plot (const std::shared_ptr<dcp::TimeDependentFunction> function,
               std::string title = "",
               std::string mode = "pause");
}

#endif
