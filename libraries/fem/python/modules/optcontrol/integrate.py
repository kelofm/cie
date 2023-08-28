from .numeric import basisfunctions
from .numeric import integration
from .numeric import functors
from .finiteelements import assembly

import math
from autograd import numpy as np




def integrateStiffnessMassRHS(n, p, length, stationarySource, k, c, basisFunctions, locationMaps, gaussPoints ):

  integrationOrder = p + 5
  numberOfDofs = p * n + 1
  K = np.zeros((numberOfDofs, numberOfDofs))
  F = np.zeros(numberOfDofs)
  M = np.zeros((numberOfDofs, numberOfDofs))
  h = length / n
  J = h / 2.0
  invJ = 2.0 / h
  
  for iElement in range( 0, n ):
    Ke = np.zeros(( p + 1, p + 1 ))
    Me = np.zeros(( p + 1, p + 1 ))
    Fe = np.zeros(p + 1)

    for iIntegrationPoint in range( 0, integrationOrder ):
      r = gaussPoints[0][iIntegrationPoint]
      w = gaussPoints[1][iIntegrationPoint]
      N = basisFunctions( r, p, 0 )
      dN = basisFunctions( r, p, 1 )
      dN *= invJ
      x = iElement * h + ( r + 1.0 ) / 2.0 * h
      f = stationarySource( x )
    
      Ke += np.transpose( dN ) * dN * k * w
      Me += np.transpose( N ) * N * c * w
      Fe += np.squeeze( N * ( f * w ) )

    Ke *= J
    Fe *= J
    Me *= J
    
    assembly.scatterMatrixInto( K, Ke, locationMaps[iElement] )
    assembly.scatterMatrixInto( M, Me, locationMaps[iElement] )
    assembly.scatterVectorInto( F, Fe, locationMaps[iElement] )

  return K, M, F

def integrateSquaredL2Norm(n, p, length, basisFunctions, locationMaps, gaussPoints, u, nTime):
  nDofs = p*n + 1
  integrationOrder = p + 5
  value = 0.0
  for iElement in range( 0, n ):
    lm = locationMaps[iElement]
    for iIntegrationPoint in range( 0, integrationOrder ):
      r = gaussPoints[0][iIntegrationPoint]
      w = gaussPoints[1][iIntegrationPoint]
      N = basisFunctions( r, p, 0 )
      for t in range(nTime):
        dofs = u[nDofs*t:nDofs*(t+1)][lm]
        value += (np.dot( N, dofs ) * w)**2
  return value


def integrateSqL2stat(n, p, length, basisFunctions, locationMaps, gaussPoints, u):
  nDofs = p*n + 1
  integrationOrder = p + 5
  value = 0.0
  for iElement in range( 0, n ):
    lm = locationMaps[iElement]
    for iIntegrationPoint in range( 0, integrationOrder ):
      r = gaussPoints[0][iIntegrationPoint]
      w = gaussPoints[1][iIntegrationPoint]
      N = basisFunctions( r, p, 0 )
      dofs = u[lm]
      value += (np.dot( N, dofs ) * w)**2
  return value


def integrateTransientRHS(n, p, length, transientSource, basisFunctions, locationMaps, gaussPoints, nTime ):
  nDofs = p*n + 1
  integrationOrder = p + 5
  h = length / n
  J = h / 2.0
  transientRHS = np.zeros(nTime*nDofs)
  for t in range(nTime):
    F = np.zeros(nDofs)
    for iElement in range( 0, n ):
      Fe = np.zeros(p + 1)
      for iIntegrationPoint in range( 0, integrationOrder ):
        r = gaussPoints[0][iIntegrationPoint]
        w = gaussPoints[1][iIntegrationPoint]
        N = basisFunctions( r, p, 0 )
        x = iElement * h + ( r + 1.0 ) / 2.0 * h
        f = transientSource( t, x )
        Fe += np.squeeze( N * ( f * w ) )
        
      Fe *= J
      assembly.scatterVectorInto( F, Fe, locationMaps[iElement] )
    transientRHS[nDofs*t:nDofs*(t+1)] = np.copy(F)
  return transientRHS  