#!/usr/bin/python

import numpy
import matplotlib.pyplot as plt

from numeric import basisfunctions

def plotUnform1DSolution( u, locationMaps, numberOfPoints, diffOrder, analyticSolutionFunction, length = 1.0, evaluateBasis = basisfunctions.IntegratedLegendre( ) ):

  n = locationMaps.shape[0]
  p = locationMaps.shape[1] - 1
  h = length / n

  uEx = [];
  uNum = [];

  N = numpy.zeros( ( numberOfPoints, p + 1 ) )
  localCoordinates = numpy.zeros( numberOfPoints )

  plotExact = analyticSolutionFunction != []

  for iResultPoint in xrange(0, numberOfPoints):

    r = -1.0 + 2.0 * iResultPoint / ( numberOfPoints - 1.0 )
    localCoordinates[iResultPoint] = r

    N[iResultPoint,:] = evaluateBasis( r, p, diffOrder )

  xVal = numpy.zeros(numberOfPoints * n)
  uNum = numpy.zeros(numberOfPoints * n)
  uEx = numpy.zeros(numberOfPoints * n)

  for iElement in xrange(0, n):      
    lm = locationMaps[iElement]
    dofs = u[lm]

    x0 = iElement * h

    for iResultPoint in xrange(0, numberOfPoints):
      r = localCoordinates[iResultPoint]
      x = x0 + ( r + 1.0 ) / 2.0 * h


      if plotExact:
        analyticSolution = analyticSolutionFunction( x )
        uEx[iElement * numberOfPoints + iResultPoint] = analyticSolution

      numericalSolution = numpy.dot( N[iResultPoint,:], dofs )
      # print("x = " + str(x) + ", numericalSolution = " + str(numericalSolution) )

      uNum[iElement * numberOfPoints + iResultPoint] = numericalSolution
      xVal[iElement * numberOfPoints + iResultPoint] = x

  if plotExact:
    p1, = plt.plot( xVal, uEx, "r", label="exact" )
    p2, = plt.plot( xVal, uNum, "b", label="numerical" )
    plt.legend(handles=[p1, p2], loc=2)      

    plt.show()

  else:
    p, = plt.plot( xVal, uNum, "b", label="numerical" )
    plt.legend(handles=[p], loc=2)      

    plt.show()

  print("_____________________")




