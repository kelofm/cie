#!/usr/bin/python

import autograd.numpy as np
from autograd import value_and_grad

from .numeric import basisfunctions
from .numeric import integration
from .numeric import functors
from .finiteelements import assembly

import math
#import numpy as np
#import matplotlib.pyplot as plt
import plotter
import integrate
import time 
import scipy.optimize

class parabolicProblem( object ):
  
  def __init__(self):    
    # ----------------- Parameters -----------------------
    
    self.l = 1.0e-2
    
    def bcfun(t):
      #if t > 6 and t < 360:
        #return 1.e4 * np.sin(1.0*(t-6)/354*np.pi)
      ##if t < 4:
        ##return 5000.0
      #else:
        #return 0.0
      return 1.0#*(0.0 + np.sin(2.0*t/80*np.pi))**2
    
    # Time
    self.nTime = 80 #time steps
    self.dt = 1.0/self.nTime   # step length
    self.theta = 0.5  # parameter for crank nicolson / theta scheme time integration
    self.pausetime = 0.005

    # Physical problem 
    # let domain be a bar of square cross section with 10mm edge length
    # norm units to meters
    self.k = 1.0#0.07#16.6 # [W / m K]
    self.c = 1.0#0.7#109.0*7874000 # [J / g K]*[g / m3]
    self.stationaryDistributedSource = lambda x: 0.0
    self.transientDistributedSource = lambda t, x: 0.0
    self.leftNeumannBoundaryFlux = bcfun #lambda t: 1.0*(0.0 + 1e2 *np.sin(1.0*t/self.nTime*np.pi))
      
    self.rightNeumannBoundaryFlux = lambda t: 0.0
    self.leftDirichletValue = 0.0
    self.rightDirichletValue = 0.0
    self.leftBC = "neumann"
    self.rightBC = "dirichlet"

    # Domain
    self.length = 1.0

    # Discretization
    self.p = 3
    self.n = 50
    self.h = self.length/self.n
    self.integrationOrder = self.p + 5
    self.nDofs = self.p * self.n + 1
    self.gaussPoints = integration.getGaussLegendreIntegrationPoints( self.integrationOrder )
    self.basisFunctions = basisfunctions.IntegratedLegendre( )

    # Boundary conditions
    self.leftDirichletBoundaryValue = 0.0
    self.rightDirichletBoundaryValue = 0.0
    self.penaltyValue = 1e8;

    # Postprocessing
    self.numberOfResultPointsPerElement = 3 + self.p * 2
    self.axisBounds = [0, self.length, -0.02, 0.02]
    self.labels = ["forw. sol.", "ref. sol.", "dual sol.", "homog. adj."]#, "app. rhs", "ref. rhs"]
    #noplot self.graphic = plotter.plotter(self.labels, self.axisBounds)
    self.plotTime = 1

    # Location Maps
    self.locationMaps = np.zeros((self.n, self.p+1)).astype(int);
    for iElement in xrange( 0, self.n ):
      locationMap = np.array([iElement, iElement + 1])
      locationMap = np.concatenate( [locationMap, np.arange( self.n + 1 + iElement * (self.p - 1), self.n + 1 + ( iElement + 1 ) * (self.p - 1) )] )
      self.locationMaps[iElement,:] = locationMap;    

    # Mass and Stiffness Matrices
    self.K, self.M, self.F = integrate.integrateStiffnessMassRHS(self.n, self.p, self.length, self.stationaryDistributedSource, self.k, self.c, self.basisFunctions, self.locationMaps, self.gaussPoints )

    # allocate arrays for concatenated solutions of all all time steps
    self.referenceSols = np.zeros(self.nDofs*self.nTime)
    self.referenceRHS = np.zeros(self.nDofs*self.nTime)
    self.forwardSols = np.zeros(self.nDofs*self.nTime)
    
    # integrate distributed right hand side
    self.referenceRHS = integrate.integrateTransientRHS(self.n, self.p, self.length, self.transientDistributedSource, self.basisFunctions, self.locationMaps, self.gaussPoints, self.nTime )

    # Boundary Conditions
    self.idLeft = 0
    self.idRight = self.n
    
    if self.rightBC == "dirichlet":
      self.K[self.idRight,self.idRight] += self.penaltyValue
      self.referenceRHS[self.idRight::self.nDofs] += self.penaltyValue * self.rightDirichletValue * np.ones(self.nTime)
    if self.leftBC == "dirichlet":
      self.K[self.idLeft,self.idLeft] += self.penaltyValue
      self.referenceRHS[self.idLeft::self.nDofs] += self.penaltyValue * self.leftDirichletValue * np.ones(self.nTime)
    if self.rightBC == "neumann":
      for t in xrange(self.nTime):
        self.referenceRHS[self.idRight+t*self.nDofs] += self.rightNeumannBoundaryFlux(t)
    if self.leftBC == "neumann":
      for t in xrange(self.nTime):
        self.referenceRHS[self.idLeft+t*self.nDofs] += self.leftNeumannBoundaryFlux(t)
        
    # Time Integration
    # precompute matrices used repeatedly for Crank-Nicolson/theta scheme integration
    self.forwLeftMatrix = (1./self.dt)*self.M + self.theta*self.K
    self.forwRightMatrix = (1./self.dt)*self.M - (1. - self.theta)*self.K
    
    # summy matrix to expand parameter vector to pde rhs
    self.dummy = np.zeros((self.nDofs*self.nTime, self.nTime))
    for i in xrange(self.nTime):
      self.dummy[self.nDofs*i, i] = 1.0

    # Compute reference solution with referenceRHS
    self.refY = np.zeros(self.nDofs) # define initial state
    t = 0 # copy initial state to solution array
    self.referenceSols[self.nDofs*t:self.nDofs*(t+1)] = self.refY[:]
    for t in xrange(1, self.nTime):
      # reference solution
      rhscn = np.dot( self.forwRightMatrix, self.refY ) + self.theta*self.referenceRHS[self.nDofs*t:self.nDofs*(t+1)] + (1-self.theta)*self.referenceRHS[self.nDofs*(t-1):self.nDofs*t]
      self.refY = np.linalg.solve(self.forwLeftMatrix, rhscn)
      self.referenceSols[self.nDofs*t:self.nDofs*(t+1)] = self.refY[:] # save solutions to long array
      
    
    
  def __call__(self, x): # ===============================================================
    
    self.approxRHS = np.dot(self.dummy, x)
    functionalValue = 0.0
    self.diffT = 0.0
    # Solve Forward Problem ---------------------------------------------------------------
    t = 0 # set initial state and copy to solution array
    y = np.zeros(self.nDofs)  # def init state
    #self.forwardSols[self.nDofs*(t):self.nDofs*(t+1)] = y[:] unnecessary
    
    #self.graphic.plot(self.n, self.p, self.h, self.basisFunctions, y, self.locationMaps, 0)
    #self.graphic.plot(self.n, self.p, self.h, self.basisFunctions, y, self.locationMaps, 1) # same init state
    #self.graphic.update('forward time: {0} of {1}'.format(t, self.nTime), self.pausetime)
      
    for t in xrange(1, self.nTime):
      rhscn = np.dot( self.forwRightMatrix, y ) + self.theta*self.approxRHS[self.nDofs*t:self.nDofs*(t+1)] + (1-self.theta)*self.approxRHS[self.nDofs*(t-1):self.nDofs*t]
      y = np.linalg.solve(self.forwLeftMatrix, rhscn)
      #self.forwardSols[self.nDofs*t:self.nDofs*(t+1)] = np.copy(y) # copy is evil
      #functionalValue += 0.5*np.sum((y-self.referenceSols[t*self.nDofs:(t+1)*self.nDofs])**2)
      functionalValue = functionalValue + 0.5*self.dt*integrate.integrateSqL2stat(self.n, self.p, self.length, self.basisFunctions, self.locationMaps, self.gaussPoints, (y-self.referenceSols[t*self.nDofs:(t+1)*self.nDofs])) + 0.5*self.dt*self.l*integrate.integrateSqL2stat(self.n, self.p, self.length, self.basisFunctions, self.locationMaps, self.gaussPoints, self.approxRHS[t*self.nDofs:(t+1)*self.nDofs])
      self.diffT = self.diffT + .5*self.dt*integrate.integrateSqL2stat(self.n, self.p, self.length, self.basisFunctions, self.locationMaps, self.gaussPoints, (y-self.referenceSols[t*self.nDofs:(t+1)*self.nDofs]))
      #if t%self.plotTime == 0:
        #self.graphic.plot(self.n, self.p, self.h, self.basisFunctions, y, self.locationMaps, 0)
        #self.graphic.plot(self.n, self.p, self.h, self.basisFunctions, self.referenceSols[self.nDofs*t:self.nDofs*(t+1)], self.locationMaps, 1)
        ##self.graphic.plot(self.n, self.p, self.h, basisFunctions, approxRHS[self.nDofs*t:self.nDofs*(t+1)], self.locationMaps, 3)
        ##self.graphic.plot(self.n, self.p, self.h, basisFunctions, self.referenceRHS[self.nDofs*t:self.nDofs*(t+1)], self.locationMaps, 4)
        #self.graphic.update('forward time: {0} of {1}'.format(t, self.nTime), self.pausetime)
    
    # Considering J() = 1/2||self.forwardSols - referenceSols||_L2 ^2
    
    return functionalValue

  def plotSim(self, x): # ===============================================================
    
    dummy = np.zeros((self.nDofs*self.nTime, self.nTime))
    for i in xrange(self.nTime):
      dummy[self.nDofs*i, i] = 1.0
      
    self.approxRHS = np.dot(dummy, x)
    
    # Solve Forward Problem ---------------------------------------------------------------
    t = 0 # set initial state and copy to solution array
    y = np.zeros(self.nDofs)  # def init state
    #self.forwardSols[self.nDofs*(t):self.nDofs*(t+1)] = y[:] unnecessary
    
    self.graphic.plot(self.n, self.p, self.h, self.basisFunctions, y, self.locationMaps, 0)
    self.graphic.plot(self.n, self.p, self.h, self.basisFunctions, y, self.locationMaps, 1) # same init state
    self.graphic.update('forward time: {0} of {1}'.format(t, self.nTime), self.pausetime)
      
    for t in xrange(1, self.nTime):
      rhscn = np.dot( self.forwRightMatrix, y ) + self.theta*self.approxRHS[self.nDofs*t:self.nDofs*(t+1)] + (1-self.theta)*self.approxRHS[self.nDofs*(t-1):self.nDofs*t]
      y = np.linalg.solve(self.forwLeftMatrix, rhscn)
      self.forwardSols[self.nDofs*t:self.nDofs*(t+1)] = np.copy(y) # copy is evil
      if t%self.plotTime == 0:
        self.graphic.plot(self.n, self.p, self.h, self.basisFunctions, y, self.locationMaps, 0)
        self.graphic.plot(self.n, self.p, self.h, self.basisFunctions, self.referenceSols[self.nDofs*t:self.nDofs*(t+1)], self.locationMaps, 1)
        self.graphic.update('forward time: {0} of {1}'.format(t, self.nTime), self.pausetime)
        
        
measurements = []

for lExp in range(-4, 5) :
  myProblem = parabolicProblem()
  myProblem.l = 10**lExp
  print "Regularization is {0}".format(myProblem.l)
  x = 0.0*np.copy(myProblem.referenceRHS[0::myProblem.nDofs])#np.zeros(myProblem.nTime)

  myProWithGrad = value_and_grad(myProblem)

  #print "value_and_grad: " + str(myProWithGrad(x))
  #print x0
  #x0[0::myProblem.nDofs] = np.ones(myProblem.nTime)*0.001

  x,f,d = scipy.optimize.fmin_l_bfgs_b(myProWithGrad, x, fprime=None, args=(), approx_grad=0, bounds=None, m=10, factr=1e0, pgtol=1e-20, iprint=1, maxfun=100, maxiter=50, disp=None, callback=None)

  #noplot myProblem.plotSim(x)
  #noplot myProblem.graphic.close()
  
  #result = scipy.optimize.fmin_cg(myProblem, x0, fprime=myProblem.getGradient, gtol=1-8, maxiter=20, full_output=True)

  #print "f: " + str(f)
  finDiffRHS = np.sum((myProblem.referenceRHS - myProblem.approxRHS)**2)
  print 'summed squared error in final rhs entries is {0}.'.format(finDiffRHS)
  finDiffT = myProblem.diffT
  print 'final l2 norm of temperature difference is {0}'.format(finDiffT)
  #print d.__str__()
  measurements.append(["{:2d}".format(lExp), "{0}".format(f), "{0}".format(finDiffT), "{0}".format(finDiffRHS)])

for m in measurements:
  print str(m)