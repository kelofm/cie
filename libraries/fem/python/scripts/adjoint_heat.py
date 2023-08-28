#!/usr/bin/python

# --- CiE Imports ---
from cie.optcontrol.numeric import basisfunctions
from cie.optcontrol.numeric import integration
import cie.optcontrol.integrate as integrate

# --- External Imports ---
import numpy as np

def call(regValue):

  # ----------------- Parameters -----------------------

  # Time
  nTime = 80   # time steps
  dt = 1.0/nTime    # step length
  theta = 0.5  # parameter for crank nicolson / theta scheme time integration

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
  l = regValue#1.0e-2   # weight for minimizing parameters as well
  maxOptIter = 50
  lowBound = 0.0
  upBound = 0.04
  functionalValue = 0.0

  # Boundary conditions
  leftDirichletBoundaryValue = 0.0
  rightDirichletBoundaryValue = 0.0
  penaltyValue = 1e8

  # Postprocessing
  numberOfResultPointsPerElement = 3 + p * 2
  axisBounds = [0, 1, -1, 1]
  labels = ["forw. sol.", "ref. sol.", "dual sol."]#, "app. rhs", "ref. rhs"]
  #noplot graphic = plotter.plotter(labels, axisBounds)
  plotTime = 100

  # Location Maps
  locationMaps = np.zeros((n, p+1)).astype(int)
  for iElement in range( 0, n ):
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
    for t in range(nTime):
      referenceRHS[idRight+t*nDofs] += rightNeumannBoundaryFlux(dt*t)
  if leftBC == "neumann":
    for t in range(nTime):
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
  for t in range(1, nTime):
    # reference solution
    rhscn = np.dot( forwRightMatrix, refY ) + theta*referenceRHS[nDofs*t:nDofs*(t+1)] + (1-theta)*referenceRHS[nDofs*(t-1):nDofs*t]
    refY = np.linalg.solve(forwLeftMatrix, rhscn)
    referenceSols[nDofs*t:nDofs*(t+1)] = refY[:] # save solutions to long array
    #noplot if t%plotTime == 0:
    #noplot     graphic.plot(n, p, h, basisFunctions, refY, locationMaps, 1)
    #noplot     graphic.update('reference time: {0} of {1}'.format(t, nTime), pausetime)
    # initialization of forward solutions
    rhscn = np.dot( forwRightMatrix, approxY ) + theta*approxRHS[nDofs*t:nDofs*(t+1)] + (1-theta)*approxRHS[nDofs*(t-1):nDofs*t]
    approxY = np.linalg.solve(forwLeftMatrix, rhscn)
    forwardSols[nDofs*t:nDofs*(t+1)] = approxY[:] # save solutions to long array


  # Considering J() = 1/2||forwardSols - referenceSols||_L2 ^2
  # possible term for smaller coefficients:  + 1/2*l*||Fapprox||_L2 ^2


  print( 'summed squared error in initial rhs entries is {0}.'.format(np.sum((referenceRHS - approxRHS)**2)) )


  # ------------------------------------------------------------------------------------------
  # Optimization Loop
  for iOpt in range(maxOptIter):
    functionalValue = 0.5*dt*integrate.integrateSquaredL2Norm(n, p, length, basisFunctions, locationMaps, gaussPoints, (forwardSols-referenceSols), nTime) + 0.5*l*dt*integrate.integrateSquaredL2Norm(n, p, length, basisFunctions, locationMaps, gaussPoints, approxRHS, nTime)
    print( 'iOpt: {0} functional value is {1}'.format(iOpt, functionalValue) )
    diffT = 0.5*dt*integrate.integrateSquaredL2Norm(n, p, length, basisFunctions, locationMaps, gaussPoints, (forwardSols-referenceSols), nTime)
    print( 'iOpt: {0} T error is {1}'.format(iOpt, diffT) )

    # construct right hand side for adjoint system
    for t in range(0,nTime):
      backwardRHS[nDofs*t:nDofs*(t+1)] = forwardSols[nDofs*(nTime-t-1):nDofs*(nTime-t)] - referenceSols[nDofs*(nTime-t-1):nDofs*(nTime-t)] # order right hand side in reverse time

    t = 0 # set initial state and copy to solution array
    dualY = np.zeros(nDofs) #forwardSols[-nDofs:] - referenceSols[-nDofs:]
    backwardSols[nDofs*t:nDofs*(t+1)] = dualY[:]

    #noplot graphic.plot(n, p, h, basisFunctions, dualY, locationMaps, 2)
    #noplot graphic.update('backward time: {0} of {1}; iOpt: {2}'.format(t, nTime, iOpt), pausetime)

    for t in range(1,nTime):   # this should run backwards?! -> ordering of entries in backSol?
      rhscn = np.dot( backwRightMatrix, dualY ) + theta*backwardRHS[nDofs*t:nDofs*(t+1)] + (1-theta)*backwardRHS[nDofs*(t-1):nDofs*t]
      dualY = np.linalg.solve(backwLeftMatrix, rhscn)
      backwardSols[nDofs*t:nDofs*(t+1)] = dualY[:]
      #noplot if t%plotTime == 0:
      #noplot   graphic.plot(n, p, h, basisFunctions, dualY, locationMaps, 2)
      #noplot   graphic.update('backward time: {0} of {1}; iOpt: {2}'.format(t, nTime, iOpt), pausetime)


    ## construct gradient by time reversing adjoint solution
    #gradient = 0.0*approxRHS[:]
    #for t in range(0,nTime):
      #gradient[nDofs*t:nDofs*(t+1)] = backwardSols[nDofs*(nTime-t-1):nDofs*(nTime-t)]
    ## only allow left boundary to change
    #for t in range(0,nTime):
      #gradient[nDofs*t+1:nDofs*(t+1)] = np.zeros(nDofs-1)

    ##print gradient.__str__()

    #approxRHS[:] = -1.0/l*gradient[:]

    projections = np.zeros(nTime)

    for t in range(0, nTime):
      value = -1.0/l*backwardSols[nDofs*(nTime-t-1)]
      if nDofs*(nTime-t-1) is 0:
        print(value)
      if value > upBound:
        value = upBound
      if value < lowBound:
        value = lowBound
      projections[t] = value

    approxRHS[::nDofs] += projections # now add the projections iteratively
    print( approxRHS[::nDofs] )

    approxY = np.zeros(nDofs)
    #noplot graphic.plot(n, p, h, basisFunctions, approxY, locationMaps, 0)
    #noplot graphic.update('forward time: {0} of {1}; iOpt: {2}'.format(0, nTime, iOpt), pausetime)
    for t in range(1, nTime):
      rhscn = np.dot( forwRightMatrix, approxY ) + theta*approxRHS[nDofs*t:nDofs*(t+1)] + (1-theta)*approxRHS[nDofs*(t-1):nDofs*t]
      approxY = np.linalg.solve(forwLeftMatrix, rhscn)
      forwardSols[nDofs*t:nDofs*(t+1)] = approxY[:] # save solutions to long array
      #noplot if t%plotTime == 0:
      #noplot   graphic.plot(n, p, h, basisFunctions, approxY, locationMaps, 0)
      #noplot   graphic.plot(n, p, h, basisFunctions, referenceSols[t*nDofs:(t+1)*nDofs], locationMaps, 1)
      #noplot   graphic.update('forward time: {0} of {1}; iOpt: {2}'.format(t, nTime, iOpt), pausetime)


    print( 'summed squared error in rhs entries is {0}.'.format(np.sum((referenceRHS - approxRHS)**2)) )



  functionalValue = 0.5*dt*integrate.integrateSquaredL2Norm(n, p, length, basisFunctions, locationMaps, gaussPoints, (forwardSols-referenceSols), nTime) + 0.5*l*dt*integrate.integrateSquaredL2Norm(n, p, length, basisFunctions, locationMaps, gaussPoints, approxRHS, nTime)
  print( 'final functional value is {0}'.format(functionalValue) )

  finDiffRHS = np.sum((referenceRHS - approxRHS)**2)

  finDiffT = 0.5*dt*integrate.integrateSquaredL2Norm(n, p, length, basisFunctions, locationMaps, gaussPoints, (forwardSols-referenceSols), nTime)
  print( 'final l2 difference in temperature is {0}'.format(finDiffT) )

  #noplot graphic.close()
  return functionalValue, finDiffT, finDiffRHS

#measurements = []
#for lExp in range(-4, 5) :
#  print( "Regularization is {0}".format(10**lExp) )
#  funVal, fDT, fDR = call(10**lExp)
#  measurements.append([funVal, fDT, fDR])
#
#for m in measurements:
#  print( str(m) )

if __name__ == "__main__":
  funVal, fDT, fDR = call( 10 )