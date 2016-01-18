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

Then, algebraic problems are tested by computing the Heaviside function of
the set defined by the points for which the solution of the linear problem
above is greater than 0.25
