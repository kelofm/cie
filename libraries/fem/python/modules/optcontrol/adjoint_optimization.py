#!/usr/bin/python

from .numeric import basisfunctions
from .numeric import integration
from .numeric import functors
from .finiteelements import assembly

import math
import numpy as np
#import matplotlib.pyplot as plt
import plotter
import integrate
import time 
import scipy.optimize as so

# ----------------- Parameters -----------------------

# Time
nTime = 80   # time steps
dt = 1.0/nTime    # step length
theta = 0.5  # parameter for crank nicolson / theta scheme time integration
pausetime = 0.005#35

# Physical problem 
k = 1.0#0.07
c = 1.0#0.7
stationaryDistributedSource = lambda x: 0.0
transientDistributedSource = lambda t, x: 0.0
leftNeumannBoundaryFlux = lambda t: 1.0#*(0.0 + np.sin(2.0*t/nTime*np.pi))**2
rightNeumannBoundaryFlux = lambda t: 0.0
leftDirichletValue = 0.0
rightDirichletValue = 0.0
leftBC = "neumann"
rightBC = "dirichlet"

# Domain
length = 1.0

# Discretization
p = 3
n = 50
h = length/n
integrationOrder = p + 5
nDofs = p * n + 1
gaussPoints = integration.getGaussLegendreIntegrationPoints( integrationOrder )
basisFunctions = basisfunctions.IntegratedLegendre( )

# Optimization
l = 1.0e-2   # weight for minimizing parameters as well
maxOptIter = 15
lowBound = 0.0
upBound = 0.04
functionalValue = 0.0

# Boundary conditions
leftDirichletBoundaryValue = 0.0
rightDirichletBoundaryValue = 0.0
penaltyValue = 1e8;

# Postprocessing
numberOfResultPointsPerElement = 3 + p * 2
axisBounds = [0, 1, -1, 1]
labels = ["forw. sol.", "ref. sol.", "dual sol."]#, "app. rhs", "ref. rhs"]
#noplot graphic = plotter.plotter(labels, axisBounds)
plotTime = 1

# Location Maps
locationMaps = np.zeros((n, p+1)).astype(int);
for iElement in xrange( 0, n ):
  locationMap = np.array([iElement, iElement + 1])
  locationMap = np.concatenate( [locationMap, np.arange( n + 1 + iElement * (p - 1), n + 1 + ( iElement + 1 ) * (p - 1) )] )
  locationMaps[iElement,:] = locationMap;    

# Mass and Stiffness Matrices
K, M, F = integrate.integrateStiffnessMassRHS(n, p, length, stationaryDistributedSource, k, c, basisFunctions, locationMaps, gaussPoints )

# allocate arrays for concatenated solutions of all all time steps
referenceSols = np.zeros(nDofs*nTime)
#referenceRHS = np.zeros(nDofs*nTime)
forwardSols = np.zeros(nDofs*nTime)
approxRHS = np.zeros(nDofs*nTime)
backwardSols = np.zeros(nDofs*nTime)
backwardRHS = np.zeros(nDofs*nTime)

# integrate distributed right hand side
referenceRHS = integrate.integrateTransientRHS(n, p, length, transientDistributedSource, basisFunctions, locationMaps, gaussPoints, nTime )

# Boundary Conditions
idLeft = 0
idRight = n

if rightBC == "dirichlet":
  K[idRight,idRight] += penaltyValue
  referenceRHS[idRight::nDofs] += penaltyValue * rightDirichletValue * np.ones(nTime)
  approxRHS[idRight::nDofs] += penaltyValue * rightDirichletValue * np.ones(nTime)
if leftBC == "dirichlet":
  K[idLeft,idLeft] += penaltyValue
  referenceRHS[idLeft::nDofs] += penaltyValue * leftDirichletValue * np.ones(nTime)
  approxRHS[idLeft::nDofs] += penaltyValue * leftDirichletValue * np.ones(nTime)
if rightBC == "neumann":
  for t in xrange(nTime):
    referenceRHS[idRight+t*nDofs] += rightNeumannBoundaryFlux(dt*t)
if leftBC == "neumann":
  for t in xrange(nTime):
    referenceRHS[idLeft+t*nDofs] += leftNeumannBoundaryFlux(t)

# Time Integration
# precompute matrices used repeatedly for Crank-Nicolson/theta scheme integration
forwLeftMatrix = (1./dt)*M + theta*K
forwRightMatrix = (1./dt)*M - (1. - theta)*K
backwLeftMatrix = (1./dt)*M + theta*K
backwRightMatrix = (1./dt)*M - (1.-theta)*K



# Compute reference solution with referenceRHS and initial forward solutions with approxRHS
refY = np.zeros(nDofs) # define initial state
approxY = np.zeros(nDofs)
t = 0 # copy initial state to solution array
referenceSols[nDofs*t:nDofs*(t+1)] = refY[:]
forwardSols[nDofs*t:nDofs*(t+1)] = approxY[:]
for t in xrange(1, nTime):
  # reference solution
  rhscn = np.dot( forwRightMatrix, refY ) + theta*referenceRHS[nDofs*t:nDofs*(t+1)] + (1-theta)*referenceRHS[nDofs*(t-1):nDofs*t]
  refY = np.linalg.solve(forwLeftMatrix, rhscn)
  referenceSols[nDofs*t:nDofs*(t+1)] = refY[:] # save solutions to long array
  #noplot if t%plotTime == 0:
  #noplot     graphic.plot(n, p, h, basisFunctions, refY, locationMaps, 1)
   #noplot    graphic.update('reference time: {0} of {1}'.format(t, nTime), pausetime)
  # initialization of forward solutions
  rhscn = np.dot( forwRightMatrix, approxY ) + theta*approxRHS[nDofs*t:nDofs*(t+1)] + (1-theta)*approxRHS[nDofs*(t-1):nDofs*t]
  approxY = np.linalg.solve(forwLeftMatrix, rhscn)
  forwardSols[nDofs*t:nDofs*(t+1)] = approxY[:] # save solutions to long array  
  

# Considering J() = 1/2||forwardSols - referenceSols||_L2 ^2
# possible term for smaller coefficients:  + 1/2*l*||Fapprox||_L2 ^2


print 'summed squared error in initial rhs entries is {0}.'.format(np.sum((referenceRHS - approxRHS)**2))


# ------------------------------------------------------------------------------------------
# Optimization Loop 
def call(x):
  approxRHS[::nDofs] = x
  
  approxY = np.zeros(nDofs)
  #noplot graphic.plot(n, p, h, basisFunctions, approxY, locationMaps, 0)
  #noplot graphic.update('forward time: {0} of {1}'.format(0, nTime), pausetime)
  for t in xrange(1, nTime):
    rhscn = np.dot( forwRightMatrix, approxY ) + theta*approxRHS[nDofs*t:nDofs*(t+1)] + (1-theta)*approxRHS[nDofs*(t-1):nDofs*t]
    approxY = np.linalg.solve(forwLeftMatrix, rhscn)
    forwardSols[nDofs*t:nDofs*(t+1)] = approxY[:] # save solutions to long array  
    #noplot if t%plotTime == 0:
    #noplot   graphic.plot(n, p, h, basisFunctions, approxY, locationMaps, 0)
    #noplot   graphic.plot(n, p, h, basisFunctions, referenceSols[t*nDofs:(t+1)*nDofs], locationMaps, 1)
    #noplot   graphic.update('forward time: {0} of {1}'.format(t, nTime), pausetime)
      
  functionalValue = 0.5*dt*integrate.integrateSquaredL2Norm(n, p, length, basisFunctions, locationMaps, gaussPoints, (forwardSols-referenceSols), nTime) + 0.5*l*dt*integrate.integrateSquaredL2Norm(n, p, length, basisFunctions, locationMaps, gaussPoints, approxRHS, nTime)
  print 'functional value is {0}'.format(functionalValue)
  
  # construct right hand side for adjoint system
  for t in xrange(0,nTime):
    backwardRHS[nDofs*t:nDofs*(t+1)] = forwardSols[nDofs*(nTime-t-1):nDofs*(nTime-t)] - referenceSols[nDofs*(nTime-t-1):nDofs*(nTime-t)] # order right hand side in reverse time
  
  t = 0 # set initial state and copy to solution array
  dualY = np.zeros(nDofs) #forwardSols[-nDofs:] - referenceSols[-nDofs:]
  #backwardSols[nDofs*t:nDofs*(t+1)] = dualY[:]
  
  #noplot graphic.plot(n, p, h, basisFunctions, dualY, locationMaps, 2)
  #noplot graphic.update('backward time: {0} of {1}'.format(t, nTime), pausetime)
  
  for t in xrange(1,nTime):   # this should run backwards?! -> ordering of entries in backSol?
    rhscn = np.dot( backwRightMatrix, dualY ) + theta*backwardRHS[nDofs*t:nDofs*(t+1)] + (1-theta)*backwardRHS[nDofs*(t-1):nDofs*t]
    dualY = np.linalg.solve(backwLeftMatrix, rhscn)
    backwardSols[nDofs*t:nDofs*(t+1)] = dualY[:]
    #noplot if t%plotTime == 0:
    #noplot   graphic.plot(n, p, h, basisFunctions, dualY, locationMaps, 2)
    #noplot   graphic.update('backward time: {0} of {1}'.format(t, nTime), pausetime)
  
  gradient = backwardSols[-1::-nDofs] * (-1.0/l)
  #gradient = gradient[::-1]
  
  return functionalValue, gradient

measurements = []

for lExp in range(-4, 5):

  l = 10**lExp
  print "Regularization is {0}".format(l)
  optRHS = so.fmin_l_bfgs_b(call, np.zeros(nTime), fprime=None, approx_grad=False, bounds=[(0, 1) for i in range(nTime)], maxiter=10, maxfun=40, disp=1, pgtol=1e-20, factr=1e0)
  
  finDiffRHS = np.sum((referenceRHS - approxRHS)**2)
  print 'summed squared error in final rhs entries is {0}.'.format(finDiffRHS)

  functionalValue = 0.5*dt*integrate.integrateSquaredL2Norm(n, p, length, basisFunctions, locationMaps, gaussPoints, (forwardSols-referenceSols), nTime) + 0.5*l*dt*integrate.integrateSquaredL2Norm(n, p, length, basisFunctions, locationMaps, gaussPoints, approxRHS, nTime)
  print 'final functional value is {0}'.format(functionalValue)

  finDiffT = 0.5*dt*integrate.integrateSquaredL2Norm(n, p, length, basisFunctions, locationMaps, gaussPoints, (forwardSols-referenceSols), nTime)
  print 'final l2 temperature difference is {0}'.format(finDiffT)
  
  measurements.append([lExp, functionalValue, finDiffT, finDiffRHS])
  
for m in measurements:
  print str(m)