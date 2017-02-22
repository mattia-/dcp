#!/usr/bin/env python
from ufl import *
set_level(DEBUG)
# define mesh dimension
cell = triangle

# define function spaces
#V = FiniteElement ("Lagrange", cell, 2)
V = VectorElement ("Lagrange", cell, 1)
Q = FiniteElement ("Lagrange", cell, 1)
T = V * Q
W = VectorElement ("Lagrange", cell, 1)
Vf = VectorElement ("Lagrange", cell, 2)
S1 = FiniteElement ("Lagrange", cell, 1)
S2 = FiniteElement ("Lagrange", cell, 2)

# define test and trial functions
#v = TestFunction (V)
vq = TestFunction (T)
v, q = split (vq)
s1 = TestFunction (S1)
s2 = TestFunction (S2)

x = SpatialCoordinate (cell)
n = FacetNormal(cell)
h_dx = Circumradius(cell)
Id = Identity(2)

nu = Constant (cell)
dt = Constant (cell)
# f = Coefficient (Vf)
f = VectorConstant (cell)
stressBelow = Coefficient (V)
beta = Constant (cell)
# wallVelocity = Coefficient (Vf)
wallVelocity = VectorConstant (cell)
gamma = Constant (cell)
cosThetaS = Constant (cell)
stabPress = Constant (cell)
stabSGCL = Constant (cell)
uno = Constant (cell)

w = Coefficient (W)
#u = Coefficient (V)
#u_old = Coefficient (V)
sol = Coefficient (T)
u,p = split (sol)
sol_old = Coefficient (T)
u_old,p_old = split (sol_old)

uDir = Coefficient (V)
#wDir = Coefficient (V)

####### - lhs

kinEn = 0.5 / dt * inner (u, u) * x[0] * dx(None, degree=5)

grav = - 1.0 / dt * inner (f, x) * x[0] * dx(None, degree=5)

gammaLength = gamma / dt * x[0] * ds(1, degree=5)
sigmaLength = - gamma * cosThetaS / dt * x[0] * ds(3, degree=5)
# sigmaLength = - 2*gamma * cosThetaS / dt * ds(3, degree=5)
sigmaLength2 = uno * x[1] * ds (3,degree=5)

kinEnFlux = 0.5 * inner (u,u) * inner (u_old,n) * x[0] * ds(2, degree=5)

gravFlux = - inner (f, x) * inner (u, n) * x[0] * ds(2,degree=5)

inflowStress = inner (stressBelow, u) * x[0] * ds(2, degree=5)

####### rhs

viscPow = 2.0 * nu * inner (sym (grad (u)), sym (grad (u))) * x[0] * dx(None, degree=5)
#viscPow = nu * inner (grad(u), grad(u)) * dx(None, degree=5)

nbc = beta*h * inner (u - wallVelocity, u) * x[0] * ds(3, degree=5)

discr = 0.5 / dt * inner (u - u_old, u - u_old) * x[0] * dx(None, degree=5)
#discrDiv = 0.5 * inner (u - u_old, u - u_old) * div(w/dt) * dx(None, degree=5)
  # w is the displacement, in this code
discrDiv = 0.5 * inner (u - u_old, u - u_old) * (Dx(x[0]*w[0],0) + x[0]*Dx(w[1],1)) / dt * dx(None, degree=5)
  # w is the displacement, in this code

####### epsilon

epsg = - 0.5 * dt * inner (f, n) * inner (w/dt, w/dt) * x[0] * ds(None, degree=5)
  # w is the displacement, in this code

tgDivw = gamma * (Dx(x[0]*(w[0]/dt),0) + x[0]*Dx(w[1]/dt,1)) * ds(1, degree=5) \
        - gamma * inner (grad(w/dt), outer(n, n)) * x[0] * ds(1, degree=5)
tgDivSol = gamma * (Dx(x[0]*(u[0]),0) + x[0]*Dx(u[1],1)) * ds(1, degree=5) \
          - gamma * inner (grad(u), outer(n, n)) * x[0] * ds(1, degree=5)

SGCLstab = stabSGCL/dt * 0.5 * gamma * ( n[0]*Dx(w[1],1) - n[1]*Dx(w[1],0) ) * ( n[0]*Dx(w[1],1) - n[1]*Dx(w[1],0) ) * x[0] * ds(1, degree=5)
SGCLstabOrig = stabSGCL * 0.5 * gamma * ( n[0]*Dx(w[1],1) - n[1]*Dx(w[1],0) ) * ( n[0]*Dx(u[0]*n[0]/n[1]+u[1],1) - n[1]*Dx(u[0]*n[0]/n[1]+u[1],0) ) * x[0] * ds(1)

bdGravw = - inner (f, x) * inner (w/dt, n) * x[0] * ds(1, degree=5)
bdGravSol = - inner (f, x) * inner (u, n) * x[0] * ds(1, degree=5)

#phiGrav = - inner (f, x) * inner (u-w/dt, n) * x[0] * ds(1, degree=5)
# NO: u^{n+1}, w^{n+1}
#phiTgDiv = gamma * (Dx(x[0]*(u[0]-w[0]/dt),0) + x[0]*Dx(u[1]-w[1]/dt,1)) * ds(1, degree=5) \
#          - gamma * inner (grad(u-w/dt), outer(n, n)) * x[0] * ds(1, degree=5)
# NO: u^{n+1}, w^{n+1}

supg = stabPress * h_dx*h_dx * inner (grad(p), grad(p)) * dx(None, degree=5)
  # NB: "*x[0]" here is neglected

####### other terms

gravPow = inner (f, u) * x[0] * dx(None, degree=5)

gDivSol = - inner (f, x) * (Dx(x[0]*u[0],0)+x[0]*Dx(u[1],1)) * dx(None, degree=5)
# gDivSol = - inner (f, x) * (Dx(x[0]*u[0],0)+x[0]*Dx(u[1],1)) * dx

tgDivSigmaw = - gamma * cosThetaS * (Dx(x[0]*(w[0]/dt),0) + x[0]*Dx(w[1]/dt,1)) * ds(3, degree=5) \
              + gamma * cosThetaS * inner (grad(w/dt), outer(n, n)) * x[0] * ds(3, degree=5)

######## GCL
intGCL1 = s1 * x[0] * dx(None, degree=5)
divGCL1 = (Dx(x[0]*w[0],0) + x[0]*Dx(w[1],1)) * s1 * dx(None, degree=5)
  # w is the displacement, in this code
intGCL2 = s2 * x[0] * dx(None, degree=5)
divGCL2 = (Dx(x[0]*w[0],0) + x[0]*Dx(w[1],1)) * s2 * dx(None, degree=5)
  # w is the displacement, in this code

forms = [kinEn, grav, gammaLength, sigmaLength,sigmaLength2, kinEnFlux, gravFlux, inflowStress, viscPow, nbc, discr, discrDiv, epsg, tgDivw, tgDivSol, SGCLstab, SGCLstabOrig, bdGravw, bdGravSol, supg, gravPow, gDivSol, tgDivSigmaw, intGCL1, divGCL1, intGCL2, divGCL2]

