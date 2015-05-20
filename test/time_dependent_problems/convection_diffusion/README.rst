This collection of files sets up the case for a test time dependent
differential problem. We will use the time dependent convection-diffusion 
equation with mixed boundary conditions on the unit square::

    \frac{\partial u}{\partial t} -k \Delta u + b \nabla u = 0 \in \Omega
    s.t.
    u (x,t) = 0 \in \partial \Omega \times (0, T]
    u (x, 0) = u_0 (x) \in \Omega

with::

    k = 0.01
    b = (1, 0)
    u_0 = exp (-[(x - 1)^2 + (y - 1)^2]/0.2)


Steps to launch the code:

1) compile::

    mkdir build
    cd build
    cmake ..
    make

2) run the executable::

    ./src/main
