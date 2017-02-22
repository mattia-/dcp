#!/usr/bin/env python
from ufl import *
set_level(DEBUG)

cell = triangle

V = FiniteElement ("Lagrange", cell, 1)
V2 = FiniteElement ("Lagrange", cell, 2)
V4 = FiniteElement ("Lagrange", cell, 4)
V6 = FiniteElement ("Lagrange", cell, 6)

v = TestFunction (V)
v2 = TestFunction (V2)
v4 = TestFunction (V4)
v6 = TestFunction (V6)

u = Coefficient (V)
datum = Coefficient (V2)
coeff = Coefficient (V2)

diffP1 = coeff * inner (u-datum, v) * ds(2, {"integration_order" : 5})
diffP2 = coeff * inner (u-datum, v2) * ds(2, {"integration_order" : 5})
diffP4 = coeff * inner (u-datum, v4) * ds(2, {"integration_order" : 5})
diffP6 = coeff * inner (u-datum, v6) * ds(2, {"integration_order" : 5})

forms = [diffP1, diffP2, diffP3]
