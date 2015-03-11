#!/usr/bin/env python
from ufl import *
set_level(DEBUG)
#/* 
# *  Copyright (C) 2014, Ivan Fumagalli, ivan.fumagalli@polimi.it
# * 
# *  This file is part of the DCP library
# *   
# *   The DCP library is free software: you can redistribute it and/or modify
# *   it under the terms of the GNU General Public License as published by
# *   the Free Software Foundation, either version 3 of the License, or
# *   (at your option) any later version.
# *
# *   The DCP library is distributed in the hope that it will be useful,
# *   but WITHOUT ANY WARRANTY; without even the implied warranty of
# *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# *   GNU General Public License for more details.
# *
# *   You should have received a copy of the GNU General Public License
# *   along with the DCP library.  If not, see <http://www.gnu.org/licenses/>. 
# */ 

# UFL file for laplace problem simulation
#
# define mesh dimension
cell = triangle

# define function spaces
V = VectorElement ("Lagrange", cell, 1)

# define test and trial functions
u = TrialFunction (V)
v = TestFunction (V)

# define coefficients
zero = Coefficient(V)

# define bilinear form
eq = inner (grad (u), grad (v)) * dx
a = eq.lhs
# define linear form 
L = eq.rhs
#L = inner(zero,v) * dx 
#inner (0*zero, v) * dx
