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

#ifndef SRC_PROBLEMS_SUBDOMAINTYPE_H_INCLUDE_GUARD
#define SRC_PROBLEMS_SUBDOMAINTYPE_H_INCLUDE_GUARD

namespace dcp
{
    /*! \brief Enumeration class used to identify domain parts
     *
     *  This class will be used to check which measure to use when setting
     *  integration subdomains on the linear form and the bilinear form.
     *  Possible values are:
     *  \li INTERNAL_CELLS: sets integration subdomains of dimension equal to that of the mesh
     *  \li BOUNDARY_FACETS: sets integration subdomains of dimension mesh_dimension - 1 on the
     *  boundary of the mesh
     *  \li INTERNAL_FACETS: sets integration subdomains of dimension mesh_dimension - 1 inside
     *  the mesh
     */
    enum class SubdomainType {INTERNAL_CELLS, BOUNDARY_FACETS, INTERNAL_FACETS};
}

#endif
