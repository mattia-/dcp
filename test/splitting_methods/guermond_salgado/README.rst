This collection of files sets up the case for a test splitting method
problem. We will use the Navier-Stokes equations with the Guermond-Salgado
splitting to study the 2D flow in the unit disk. 
See J.-L. Guermond, A. Salgado, A splitting method for incompressible flows with variable
density ..., J. Comput. Phys. (2009), doi:10.1016/j.jcp.2008.12.036, section 5.1

Steps to launch the code:

1) compile::

    mkdir build
    cd build
    cmake ..
    make

You may have to adjust the directories where cmake looks for header files and libraries if you installed 
DCP in a non-default location. To do this, use the variables INCLUDE_DIRS and LINK_DIRS when invoking cmake

2) run the executable::

    ./src/main
