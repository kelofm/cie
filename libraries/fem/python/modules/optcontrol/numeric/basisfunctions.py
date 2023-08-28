#!/usr/bin/python

import math
from autograd import numpy
import matplotlib.pyplot as plt

from .integration import *

class GridTensorProduct( object ):

  def __init__( self, basisFunctions ):
    self.basisFunctions = basisFunctions


  def evaluate( self, rst, p, diffOrder, numberOfFields = None ):

    # ----------------- Input Checking ----------------------

    # check input type of rst
    if not hasattr( rst, "__len__"):
      ValueError( 'rst has to be a list type!' )

    spaceDimension = len( rst )

    if not hasattr( p, "__len__"):
      n = numberOfFields

      if n == None:
        n = 1

      p = numpy.full( ( n, spaceDimension ), p, dtype = numpy.int )

    if not hasattr( p, "shape" ):
      p = numpy.array( p )

    if numberOfFields == None:
      numberOfFields = p.shape[0];
    
    if p.shape != ( numberOfFields, spaceDimension ):
      raise ValueError( "Order template dimensions are not consistent!" )

    if diffOrder < 0 or diffOrder > 2:
      raise ValueError( "Diff order not supported!" )

    # -------------- Compute 1D basis functions -------------------

    evaluationOf1DBasis = []
  
    for i in range( spaceDimension ):

      evaluationOf1DBasisInI = {};  
      for polynomialOrder in p[:,i]:

        if not polynomialOrder in evaluationOf1DBasisInI:
          evaluationOf1DBasisInIForP0 = self.basisFunctions( rst[i][0], polynomialOrder, 0 )

          numberOf1DBasisFunctions = evaluationOf1DBasisInIForP0.shape[1]
          numberOfEvaluationPoints = len( rst[i] );

          evaluationOf1DBasisInIForP = numpy.zeros( ( diffOrder + 1, numberOfEvaluationPoints, numberOf1DBasisFunctions ) );
          evaluationOf1DBasisInIForP[0, 0] = evaluationOf1DBasisInIForP0

          # Fill remaining diff orders for first evaluation coordinate
          for j in range( 1, diffOrder + 1 ):
            evaluationOf1DBasisInIForP[j, 0] = self.basisFunctions( rst[i][0], polynomialOrder, j )

          # Compute the other evaluation coordinates for all diff orders
          for j in range( 1, numberOfEvaluationPoints ):
            for k in range( 0, diffOrder + 1 ):
              evaluationOf1DBasisInIForP[k, j] = self.basisFunctions( rst[i][j], polynomialOrder, k )

          evaluationOf1DBasisInI[polynomialOrder] = evaluationOf1DBasisInIForP
  
      evaluationOf1DBasis.append( evaluationOf1DBasisInI )

    # -------------- Compute Tensor product -------------------

    differentPolynomialOrders = numpy.any( numpy.diff(p, axis=0) != 0 )

    shapes = [];
    
    for iField in range( numberOfFields ):    
      diffOrders = [];
      
      shapesInR = evaluationOf1DBasis[0][p[iField, 0]]
      shapesInS = evaluationOf1DBasis[1][p[iField, 1]]

      # Form 4d array with evaluation points in dimensions 0 & 1 and tensor product in 2 & 3 
      N = ( shapesInR[0, :, numpy.newaxis, :, numpy.newaxis] * 
            shapesInS[0, numpy.newaxis, :, numpy.newaxis, :] )

      # merge last two dimensions
      N.shape = ( N.shape[0], N.shape[1], N.shape[2] * N.shape[3] )

      diffOrders.append( N )

      if diffOrder > 0:

        dN = numpy.zeros( ( len(rst[0]), len(rst[1]), 2, shapesInR.shape[2], shapesInS.shape[2] ) )
        
        dN[:,:,0,:,:] = ( shapesInR[1, :, numpy.newaxis, :, numpy.newaxis] * 
                          shapesInS[0, numpy.newaxis, :, numpy.newaxis, :] )        
        dN[:,:,1,:,:] = ( shapesInR[0, :, numpy.newaxis, :, numpy.newaxis] * 
                          shapesInS[1, numpy.newaxis, :, numpy.newaxis, :] )

        # merge last two dimensions
        dN.shape = ( dN.shape[0], dN.shape[1], dN.shape[2], dN.shape[3] * dN.shape[4] )

        diffOrders.append( dN )

      if diffOrder > 1:

        ddN = numpy.zeros( ( len(rst[0]), len(rst[1]), 3, shapesInR.shape[2], shapesInS.shape[2] ) )
        
        ddN[:,:,0,:,:] = ( shapesInR[2, :, numpy.newaxis, :, numpy.newaxis] * 
                           shapesInS[0, numpy.newaxis, :, numpy.newaxis, :] )        
        ddN[:,:,1,:,:] = ( shapesInR[0, :, numpy.newaxis, :, numpy.newaxis] * 
                           shapesInS[2, numpy.newaxis, :, numpy.newaxis, :] )
        ddN[:,:,2,:,:] = ( shapesInR[1, :, numpy.newaxis, :, numpy.newaxis] * 
                           shapesInS[1, numpy.newaxis, :, numpy.newaxis, :] )        

        # merge last two dimensions
        ddN.shape = ( ddN.shape[0], ddN.shape[1], ddN.shape[2], ddN.shape[3] * ddN.shape[4] )

        diffOrders.append( ddN )

      shapes.append( diffOrders )

      # if the other field components are equal, we can copy
      if differentPolynomialOrders == False:
        for jField in range( 1, numberOfFields ):
          diffOrdersCopy = []
          for iDiffOrder in range( diffOrder + 1 ):
            diffOrdersCopy.append( numpy.copy( diffOrders[iDiffOrder] ) ) 
          shapes.append( diffOrdersCopy )
        break

    return shapes

  def createNumberingForOrderTemplate( self, p ):
    
    numberOfFields, spaceDim = p.shape

    numbering = [];

    for iFieldComponent in range( numberOfFields ):
      degreesForFieldComponent = p[0, :] + 1;

      if spaceDim == 1:
        numering.append( numpy.arange( degreesForFieldComponent[0, 0] ) )  
      
      if spaceDim == 2:
        mapping = numpy.zeros( degreesForFieldComponent, dtype=numpy.int )

        pMinusOneR = degreesForFieldComponent[0] - 2;
        pMinusOneS = degreesForFieldComponent[1] - 2;

        # nodal modes
        mapping[0:2, 0] = [0, 2]
        mapping[0:2, 1] = [1, 3]

        # edge modes
        index = 4;
        mapping[2:, 0] = numpy.arange( index             , index +     pMinusOneR )
        mapping[2:, 1] = numpy.arange( index + pMinusOneR, index + 2 * pMinusOneR )

        index += 2 * pMinusOneR;
        mapping[0, 2:] = numpy.arange( index             , index +     pMinusOneS )
        mapping[1, 2:] = numpy.arange( index + pMinusOneS, index + 2 * pMinusOneS )
        
        index += 2 * pMinusOneS;

        volumeModeIndices = numpy.ix_( numpy.arange( 2, degreesForFieldComponent[0] ),   
                                       numpy.arange( 2, degreesForFieldComponent[0] ) )
        mapping[ volumeModeIndices ] = numpy.reshape( numpy.arange( index, index + pMinusOneR * pMinusOneS ), ( pMinusOneS, pMinusOneR ) ).T

        numbering.append( mapping )

    return numbering

## Integrated legendre basis functions
class IntegratedLegendre( object ):

  ## Constructor, doing nothing, just initializing with zero.
  #  @param self               The object pointer.
  def __init__( self ):
    
    self.maxFactorIndex = 0
    
    self.factors = numpy.zeros( 0 )
    self.oneOverJ = numpy.zeros( 0 )
    self.twoJminusOne = numpy.zeros( 0 )
    self.JminusOne = numpy.zeros( 0 )

  ## Call to IntegratedLegendre functor to evaluate shape functions.
  #  @param self      The object pointer.
  #  @param r         The local 1d coordinate in [-1, 1].
  #  @param pMax      The polynomial degree up to which evaluation should be performed.
  #  @param diffOrder The order of differention.
  def __call__( self, r, pMax, differentialOrder ):

    # check if factors are there, compute if not
    if pMax > self.maxFactorIndex - 1:
      self._computeFactorsUpTo( pMax + 1 )  

    if differentialOrder == 0:

      # Julian added some lines to ease least dquares fitting, i.e. let this return values for a single linear function
      if pMax < 1:
        pMax = 1
      #if pMax == 0:
        #N = numpy.zeros((1,1))
        #N[0,0] = 0.5 * ( 1.0 - r )
        #return N
      
      # allocate memory for shapes and set linear modes
      N = numpy.zeros( ( 1, pMax + 1 ) )
      N[0, 0] = 0.5 * ( 1.0 - r )
      N[0, 1] = 0.5 * ( 1.0 + r )

      # allocate memory for the three last Legendre shape functions
      L = [ 1, r, 0.0 ]

      for j in range( 2, pMax + 1 ):

        # compute the next Legendre function
        L[2] = self.oneOverJ[j] * ( self.twoJminusOne[j] * r * L[1] - self.JminusOne[j] * L[0] )

        # compute the next integrated Legrendre shape functions
        N[0, j] = self.factors[j] * ( L[2] - L[0] )

        # and finally move the Legendre shape functions
        L[0] = L[1]
        L[1] = L[2]
 
      return N

    elif differentialOrder == 1:

      # allocate memory for shape derivatives and set linear modes
      dN = numpy.zeros( ( 1, pMax + 1 ) )
      dN[0, 0] = -0.5
      dN[0, 1] =  0.5

      # allocate memory for the three last Legendre shape functions and their derivatives
      L  = [ 1, r, 0.0 ]
      dL = [ 0, 1, 0.0 ]

      for j in range( 2, pMax + 1 ):

        # compute the next Legendre function
        L[2] = self.oneOverJ[j] * ( self.twoJminusOne[j] * r * L[1] - self.JminusOne[j] * L[0] )

        # and its derivative
        dL[2] = self.oneOverJ[j] * ( self.twoJminusOne[j] * ( L[1] + r * dL[1] ) - self.JminusOne[j] * dL[0] )

        # compute the next integrated Legrendre shape functions
        dN[0, j] = self.factors[j] * ( dL[2] - dL[0] )

        # and finally move the Legendre shape functions
        L[0] = L[1]
        L[1] = L[2]

        # and and their derivatives
        dL[0] = dL[1]
        dL[1] = dL[2]

      return dN

    elif ( differentialOrder == 2 ):

      # allocate memory for shape derivatives
      ddN = numpy.zeros( ( 1, pMax + 1 ) )

      # allocate memory for the three last Legendre shape functions and their derivatives
      L   = [ 1, r, 0.0 ]
      dL  = [ 0, 1, 0.0 ]
      ddL = [ 0, 0, 0.0 ]

      for j in range( 2, pMax + 1 ):

        # compute the next Legendre function
        L[2] = self.oneOverJ[j] * ( self.twoJminusOne[j] * r * L[1] - self.JminusOne[j] * L[0] )

        # and its first derivative
        dL[2] = self.oneOverJ[j] * ( self.twoJminusOne[j] * ( L[1] + r * dL[1] ) - self.JminusOne[j] * dL[0] )

        # and its second derivative
        ddL[2] = self.oneOverJ[j] * ( self.twoJminusOne[j] * ( 2.0 * dL[1] + r * ddL[1] ) - self.JminusOne[j] * ddL[0] )

        # compute the next integrated Legrendre shape functions
        ddN[0, j] = self.factors[j] * ( ddL[2] - ddL[0] )

        # and finally move the Legendre shape functions
        L[0] = L[1]
        L[1] = L[2]

        # and and their derivatives
        dL[0] = dL[1]
        dL[1] = dL[2]

        ddL[0] = ddL[1]
        ddL[1] = ddL[2]

      return ddN

    else:
      raise ValueError('Differential order not supported!')

  ## Internal functions to precompute values needed to evaluate the basis functions.
  #  @param self               The object pointer.
  #  @param computeFactorsUpTo The polynomial degree up to which evaluation should be performed.
  def _computeFactorsUpTo( self, index ):
    
    self.maxFactorIndex = index
    
    self.factors = numpy.zeros( index )
    self.oneOverJ = numpy.zeros( index )
    self.twoJminusOne = numpy.zeros( index )
    self.JminusOne = numpy.zeros( index )

    for j in range( 2, index ):
      self.factors[j] = 1.0 / math.sqrt( 4.0 * j - 2.0 )
      self.twoJminusOne[j] = 2.0 * j - 1.0
      self.JminusOne[j] = j - 1.0
      self.oneOverJ[j] = 1.0 / j

## Lagrange basis functions
class Lagrange( object ):

  def __init__( self, samplePointProvider = None ):
    
    if samplePointProvider == None:
      self.samplePointProvider = integration.IntegrationPointProvider( integration.getTrapezoidalIntegrationPoints )
    else:
      self.samplePointProvider = samplePointProvider

    self.coefficients = { }

  def __call__( self, r, p, differentialOrder ):

    if p < 1:
      raise ValueError( 'Polynomial degree has to be > 0!' )
    elif r < -1.0 or r > 1.0:
      raise ValueError( 'Local coordinate has to be in [-1, 1]!' )    

    if not p in self.coefficients:
      self._computeCoefficients( p )

    C = self.coefficients[p]
    N = numpy.ones( ( 1, p + 1 ) )

    if differentialOrder == 0:   

      for i in range( p + 1 ):
        N[0, i] = C[i, p]

        for j in range( p - 1, -1, -1 ):
          N[0, i] = N[0, i] * r + C[i, j]

    elif differentialOrder == 1:  

      for i in range( p + 1 ):
        N[0, i] = C[i, p] * float(p)
      
        for j in range( p - 1, 0, -1 ):
          N[0, i] = N[0, i] * r + C[i, j] * float( j )

    elif differentialOrder == 2:

      for i in range( p + 1 ):
        N[0, i] = C[i, p] * float(p) * ( float(p) - 1.0 )
      
        for j in range( p - 1, 1, -1 ):
          N[0, i] = N[0, i] * r + C[i, j] * float( j ) * ( float( j ) - 1.0 )

    else:
      raise ValueError( 'Differential order not supported!' )

    return N

  def _computeCoefficients( self, n ):
    
    mySamplingPoints = self.samplePointProvider( n + 1 )[0, :]

    tmpCoefficients = numpy.zeros( n + 1 )
    tmpCoefficients[n] = -mySamplingPoints[0]

    myCoefficients = numpy.zeros( ( n + 1, n + 1 ) )

    for i in range( 1, n + 1 ):
      for j in range( n - i, n ):
        tmpCoefficients[j] -= mySamplingPoints[i] * tmpCoefficients[j + 1]

      tmpCoefficients[n] -= mySamplingPoints[i]

    # Calculate the actual coefficients
    for index in range( n + 1 ):
      for j in range( n + 1 ):

        phi = n + 1
        for k in range( n, 0, -1 ):
          phi = k * tmpCoefficients[k] + mySamplingPoints[j] * phi

        if j == index:
          ff = 1.0 / phi

        else:
          ff = 0.0

        b = 1.0

        for k in range( n, -1, -1 ):
          myCoefficients[index][k] += b * ff
          b = tmpCoefficients[k] + mySamplingPoints[j] * b

    self.coefficients[ n ] = myCoefficients

def plotIntegratedLegendrePolynomials( pMin, pMax, derivative = 0 ):
  
  nPoints = 200

  nFunctions = pMax - pMin + ( pMin == 1 ) + 1
  minIndex = pMin - (pMin == 1)

  evaluateFunction = IntegratedLegendre()

  r = numpy.linspace(-1, 1, nPoints)

  f = numpy.zeros((nFunctions, nPoints))

  for i in range(0, nPoints):
    N = evaluateFunction(r[i], pMax, derivative)
    f[:,i] = N[:,minIndex:(minIndex+nFunctions)]
  
  for i in range(0, nFunctions):
    plt.plot(r, f[i,:])

  plt.show( )
  


def evaluateLegendrePolynomialRecursive( n, r ):
    if n > 1:
      return ( 1.0 / n * ( ( 2.0 * n - 1.0 ) * r * evaluateLegendrePolynomialRecursive( n - 1.0, r ) - 
                                     ( n - 1.0 ) * evaluateLegendrePolynomialRecursive( n - 2.0, r ) ) )

    elif n == 1:
      return ( r )

    else: # n == 0 
      return ( 1.0 )

# def evaluateLegendrePolynomialDerivativeRecursive( n, r ): 
#     if n > 1:
#       return ( ( 2.0 * n - 1.0 ) * evaluateLegendrePolynomialRecursive( n - 1.0, r )
#                                  + evaluateLegendrePolynomialDerivativeRecursive( n - 2.0, r ) )
# 
#     elif n == 1:
#       return 1.0
# 
#     else: # n == 0 
#       return 0.0
# 
# 
# def evaluateLegendrePolynomialSecondDerivativeRecursive( n, x ):
#     if n > 1:
#       return ( ( 2.0 * n - 1.0 ) * evaluateLegendrePolynomialDerivativeRecursive( n - 1.0, x )
#                                  + evaluateLegendrePolynomialSecondDerivativeRecursive( n - 2.0, x ) )
# 
#     else: # n == 0 || n == 1
#       return 0.0
# 
# def evaluateLegendrePolynomial(r, p, derivative):
# 
#   r2 = r*r
# 
#   if derivative == 0:
# 
#     if   p == 0:  return ( 1.0 )
#     elif p == 1:  return ( r )
#     elif p == 2:  return ( 3.0/2.0*r2-1.0/2.0 )
#     elif p == 3:  return ( (-3.0/2.0+5.0/2.0*r2 )*r )
#     elif p == 4:  return ( 3.0/8.0+(-15.0/4.0+35.0/8.0*r2 )*r2 )
#     elif p == 5:  return ( (15.0/8.0+(-35.0/4.0+63.0/8.0*r2 )*r2 )*r )
#     elif p == 6:  return ( -5.0/16.0+(105.0/16.0+(-315.0/16.0+231.0/16.0*r2 )*r2 )*r2 )
#     elif p == 7:  return ( (-35.0/16.0+(315.0/16.0+(-693.0/16.0+429.0/16.0*r2 )*r2 )*r2 )*r )
#     elif p == 8:  return ( 35.0/128.0+(-315.0/32.0+(3465.0/64.0+(-3003.0/32.0+6435.0/128.0
#                            *r2 )*r2 )*r2 )*r2 )
#     elif p == 9:  return ( (315.0/128.0+(-1155.0/32.0+(9009.0/64.0+(-6435.0/32.0+12155.0/128.0
#                            * r2 )*r2 )*r2 )*r2 )*r )
#     elif p == 10: return ( -63.0/256.0+(3465.0/256.0+(-15015.0/128.0+(45045.0/128.0+(
#                            -109395.0/256.0+46189.0/256.0*r2 )*r2 )*r2 )*r2 )*r2 )
#     elif p == 11: return ( (-693.0/256.0+(15015.0/256.0+(-45045.0/128.0+(109395.0/128.0
#                            +(-230945.0/256.0+88179.0/256.0*r2 )*r2 )*r2 )*r2 )*r2 )*r )
#     elif p == 12: return ( 231.0/1024.0+(-9009.0/512.0+(225225.0/1024.0+(-255255.0/256.0
#                            +(2078505.0/1024.0+(-969969.0/512.0+676039.0/1024.0*r2 )*r2 )*r2 
#                            )*r2 )*r2)*r2 )
#     elif p == 13: return ( (3003.0/1024.0+(-45045.0/512.0+(765765.0/1024.0+(-692835.0/256.0
#                            +(4849845.0/1024.0+(-2028117.0/512.0+1300075.0/1024.0*r2 )*r2)
#                            *r2 )*r2 )*r2 )*r2 )*r )
#     elif p == 14: return ( -429.0/2048.0+(45045.0/2048.0+(-765765.0/2048.0+(4849845.0/2048.0
#                            +(-14549535.0/2048.0+(22309287.0/2048.0+(-16900975.0/2048.0
#                            +5014575.0/2048.0*r2)*r2 )*r2 )*r2 )*r2 )*r2 )*r2 )
#     elif p == 15: return ( (-6435.0/2048.0+(255255.0/2048.0+(-2909907.0/2048.0+(14549535.0
#                            /2048.0+(-37182145.0/2048.0+(50702925.0/2048.0+(-35102025.0
#                            /2048.0+9694845.0/2048.0*r2 )*r2 )*r2 )*r2 )*r2 )*r2 )*r2 )*r )
#     elif p == 16: return ( 6435.0/32768.0+(-109395.0/4096.0+(4849845.0/8192.0
#                            +(-20369349.0/4096.0+(334639305.0/16384.0+(-185910725.0/4096.0
#                            +(456326325.0/8192.0+(-145422675.0/4096.0+300540195.0/32768.0*r2 )
#                            *r2 )*r2 )*r2 )*r2 )*r2 )*r2)*r2 )
#     elif p == 17: return ( (109395.0/32768.0+(-692835.0/4096.0+(20369349.0/8192.0
#                            +(-66927861.0/4096.0+(929553625.0/16384.0+(-456326325.0/4096.0
#                            +(1017958725.0/8192.0+(-300540195.0/4096.0+583401555.0/32768.0*r2 )
#                              *r2 )*r2 )*r2 )*r2 )*r2 )*r2 )*r2 )*r )
#     elif p == 18: return ( -12155.0/65536.0+(2078505.0/65536.0+(-14549535.0/16384.0
#                            +(156165009.0/16384.0+(-1673196525.0/32768.0+(5019589575.0
#                            /32768.0+(-4411154475.0/16384.0+(4508102925.0/16384.0+(
#                            -9917826435.0/65536.0+2268783825.0/65536.0*r2)*r2 )*r2 )*r2 )*r2 )
#                            *r2 )*r2 )*r2 )*r2 )
#     elif p == 19: return ( (-230945.0/65536.0+(14549535.0/65536.0+(-66927861.0/16384.0
#                            +(557732175.0/16384.0+(-5019589575.0/32768.0+(13233463425.0
#                            /32768.0+(-10518906825.0/16384.0+(9917826435.0/16384.0+(
#                            -20419054425.0/65536.0+4418157975.0/65536.0*r2 )*r2 )*r2 )*r2 )
#                            *r2 )*r2 )*r2 )*r2 )*r2 )*r )
#     elif p == 20: return ( 46189.0/262144.0+(-4849845.0/131072.0+(334639305.0/262144.0+(
#                            -557732175.0/32768.0+(15058768725.0/131072.0+(-29113619535.0
#                            /65536.0+(136745788725.0/131072.0+(-49589132175.0/32768.0
#                            +(347123925225.0/262144.0+(-83945001525.0/131072.0+34461632205.0
#                            /262144.0*r2 )*r2 )*r2 )*r2 )*r2 )*r2 )*r2 )*r2 )*r2 )*r2  )
#     else:
#       return evaluateLegendrePolynomialRecursive( p, r )
# 
#   elif derivative == 1:
# 
#     if   p ==  0: return ( 0.0 )
#     elif p ==  1: return ( 1.0 )            
#     elif p ==  2: return ( 3.0*r )                
#     elif p ==  3: return ( 15.0/2.0*r2-3.0/2.0 )        
#     elif p ==  4: return ( (-15.0/2.0+35.0/2.0*r2 )*r )        
#     elif p ==  5: return ( 15.0/8.0+(-105.0/4.0+315.0/8.0*r2 )*r2 )        
#     elif p ==  6: return ( (105.0/8.0+(-315.0/4.0+693.0/8.0*r2 )*r2 )*r )        
#     elif p ==  7: return ( -35.0/16.0+(945.0/16.0+(-3465.0/16.0+3003.0/16.0*r2 )*r2 )*r2)
#     elif p ==  8: return ( (-315.0/16.0+(3465.0/16.0+(-9009.0/16.0+6435.0/16.0*r2 )*r2 )*
#                            r2 )*r )
#     elif p ==  9: return ( 315.0/128.0+(-3465.0/32.0+(45045.0/64.0+(-45045.0/32.0+109395.0/
#                            128.0*r2 )*r2 )*r2 )*r2 )
#     elif p == 10: return ( (3465.0/128.0+(-15015.0/32.0+(135135.0/64.0+(-109395.0/32.0+
#                            230945.0/128.0*r2 )*r2 )*r2 )*r2 )*r )        
#     elif p == 11: return ( -693.0/256.0+(45045.0/256.0+(-225225.0/128.0+(765765.0/128.0+(
#                            -2078505.0/256.0+969969.0/256.0*r2 )*r2 )*r2 )*r2 )*r2 )        
#     elif p == 12: return ( (-9009.0/256.0+(225225.0/256.0+(-765765.0/128.0+(2078505.0/128.0+(
#                            -4849845.0/256.0+2028117.0/256.0*r2 )*r2 )*r2 )*r2 )*r2 )*r )
#     elif p == 13: return ( 3003.0/1024.0+(-135135.0/512.0+(3828825.0/1024.0+(-4849845.0/256.0
#                            +(43648605.0/1024.0+(-22309287.0/512.0+16900975.0/1024.0*r2 )*r2 )
#                            *r2 )*r2)*r2 )*r2 )
#     elif p == 14: return ( (45045.0/1024.0+(-765765.0/512.0+(14549535.0/1024.0+(-14549535.0/
#                            256.0+(111546435.0/1024.0+(-50702925.0/512.0+35102025.0/1024.0
#                            *r2 )*r2)*r2 )* r2 )*r2 )*r2 )*r )        
#     elif p == 15: return ( -6435.0/2048.0+(765765.0/2048.0+(-14549535.0/2048.0+(101846745.0/
#                            2048.0+(-334639305.0/2048.0+(557732175.0/2048.0+(-456326325.0/
#                            2048.0+145422675.0/2048.0*r2 )*r2 )*r2 )*r2 )*r2 )*r2 )*r2 )
#     elif p == 16: return ( (-109395.0/2048.0+(4849845.0/2048.0+(-61108047.0/2048.0+(
#                            334639305.0/2048.0+(-929553625.0/2048.0+(1368978975.0/2048.0+(-
#                            1017958725.0/2048.0+300540195.0/2048.0*r2 )*r2 )*r2 )*r2 )*r2 )*
#                            r2 )*r2 )*r )        
#     elif p == 17: return ( 109395.0/32768.0+(-2078505.0/4096.0+(101846745.0/8192.0+(-
#                            468495027.0/4096.0+(8365982625.0/16384.0+(-5019589575.0/4096.0+(
#                            13233463425.0/8192.0+(-4508102925.0/4096.0+9917826435.0/32768.0*
#                            r2 )*r2 )*r2 )*r2 )*r2)*r2 ) *r2 )*r2 )        
#     elif p == 18: return ( (2078505.0/32768.0+(-14549535.0/4096.0+(468495027.0/8192.0+(-
#                            1673196525.0/4096.0+(25097947875.0/16384.0+(-13233463425.0/
#                            4096.0+(31556720475.0/8192.0+(-9917826435.0/4096.0+20419054425.0/
#                            32768.0*r2 )*r2)*r2 )*r2 )*r2 )*r2 )*r2 )*r2 )*r )  
#     elif p == 19: return ( -230945.0/65536.0+(43648605.0/65536.0+(-334639305.0/16384.0+(
#                            3904125225.0/16384.0+(-45176306175.0/32768.0+(145568097675.0/
#                            32768.0+(-136745788725.0/16384.0+(148767396525.0/16384.0+(-
#                            347123925225.0/65536.0+83945001525.0/65536.0*r2 )*r2 )*r2 )*r2 )*
#                            r2 )*r2 )*r2 )*r2 )*r2 )        
#     elif p == 20: return ( (-4849845.0/65536.0+(334639305.0/65536.0+(-1673196525.0/16384.0+(
#                            15058768725.0/16384.0+(-145568097675.0/32768.0+(410237366175.0/
#                            32768.0+(-347123925225.0/16384.0+(347123925225.0/16384.0+(-
#                            755505013725.0/65536.0+172308161025.0/65536.0*r2 )*r2 )*r2 )*r2 )*
#                            r2 )*r2 )*r2 )*r2 )*r2 )*r )
#     else:
#      return  evaluateLegendrePolynomialDerivativeRecursive( p, r )
# 
#   elif derivative == 2:
#     return evaluateLegendrePolynomialSecondDerivativeRecursive( p, r )
# 
#   else:
#     raise ValueError('Unsupported diff order!')
# 
# 
# def evaluateIntegratedLegendrePolynomial(r, i, derivative):
# 
#   if derivative == 0:
#     if i == 0:
#       return 0.5*(1-r)
# 
#     elif i == 1:
#       return 0.5*(1+r)
# 
#     else:
#       return ( evaluateLegendrePolynomial(r,i, 0) - evaluateLegendrePolynomial(r,i-2, 0) ) / math.sqrt(4*i-2)
# 
#   elif derivative == 1:
#     if i == 0:
#       return -0.5
# 
#     elif i == 1:
#       return 0.5
# 
#     else:
#       return ( evaluateLegendrePolynomial(r,i, 1) - evaluateLegendrePolynomial(r,i-2, 1) ) / math.sqrt(4*i-2)
#   
#   elif derivative == 2:
#     if i < 2:
#       return 0.0
#     else:
#       return( evaluateLegendrePolynomial(r,i, 2) - evaluateLegendrePolynomial(r,i-2, 2) ) / math.sqrt(4*i-2)
# 
#   else:
#     raise Exeption('Unkown derivative order')
# 
# def evaluate1DIntegratedLegendreBasis( r, p, derivative ):
#   
#   N = numpy.zeros( (1, p+1) )
# 
#   for i in range(0, p+1 ):
#     N[0, i] = evaluateIntegratedLegendrePolynomial( r, i, derivative )
# 
#   return N





