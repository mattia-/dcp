This collection of files sets up the case for a test nonlinear
differential problem.
We will use the stationary Navier-Stokes equations to study the flow
around a cylinder.

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
