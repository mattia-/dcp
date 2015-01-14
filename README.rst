Author: Mattia Tamellini, 2014  <mattia.tamellini@gmail.com>

===
DCP
===

DCP is a library for optimal control problems based on dolfin/FEniCS: http://fenicsproject.org/support/

Compilation
===========

To build DCP, run::

    mkdir build
    cd build
    cmake ..
    make

various options can be passed to ``cmake``:  

- ``C_COMPILER`` and ``CXX_COMPILER`` : the C and C++ compiler to be used

- ``CMAKE_BUILD_TYPE`` : sets the compilation type, to allow for some optimization at
  compile time (see CMake documentation for more details)

- ``LIB_INSTALL_DIR`` : directory for the installation of the library. 
  Default value: the directory `lib` in the root directory of the source tree

- ``HEADER_INSTALL_DIR`` : directory for the installation of the header files. 
  Default value: the directory `include` in the root directory of the source tree 

All option can be set from command line using the flag ``-D`` : ::

    cmake .. -DVARIABLE_NAME=VARIABLE_VALUE   


CMake generated makefiles support parallel compilation through the flag ``-j``. 
For example, to compile using 4 threads use the command: ::
    
    make -j 4


Installation
============

From the build directory, launch::
    
    make install


Test
====

To run a test, go to the ``test`` directory and follow the instructions of the file README in there
  
License
=======

DCP is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

DCP is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with DCP. If not, see <http://www.gnu.org/licenses/>.
