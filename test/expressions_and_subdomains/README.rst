This collection of files sets up the case for a test of the classes
dcp::VariableExpression and dcp::Subdomain, and in particular their
use through callable objects.
We will solve a linear differential problem, namely the Poisson
equation with mixed boundary conditions on the unit square::

    -k \Delta u = f \in \Omega
    s.t.
    \nabla u \cdot n = g \in \Gamma_N
    u = h \in \partial \Omega \setminus \Gamma_N

with::

    \Gamma_N = {(x, y) \in \Omega s.t. x \in [0, 1] and y = 0}
    k = 1
    f = 1
    g = x^2
    h = x+y
