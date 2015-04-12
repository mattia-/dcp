This collection of files sets up the case for a test linear differential
problem. We will use the Poisson equation with mixed boundary conditions
on the unit square::

    -k \Delta u = f \in \Omega
    s.t.
    \nabla u \cdot n = g \in \Gamma_N
    u = h \in \partial \Omega \setminus \Gamma_N

with::

    \Gamma_N = {(x, y) \in \Omega s.t. x \in [0, 1] and y = 0}
    k = 1
    f = 1
    g = 1
    h = 0


Steps to launch the code:

1) compile::

    mkdir build
    cd build
    cmake ..
   
You may have to adjust the directories where cmake looks for header files and libraries if you installed 
DCP in a non-default location. To do this, use the variables INCLUDE_DIRS and LINK_DIRS when invoking cmake

2) run the executable::

    ./src/main
