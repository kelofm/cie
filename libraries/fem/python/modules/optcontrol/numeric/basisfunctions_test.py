#!/usr/bin/python

# python imports
import unittest
import numpy

# numeric imports
import basisfunctions

class TestLagrangeBasisFunctions(unittest.TestCase):

  def setUp( self ):
    self.tolerance = 1e-10

  def test_difforder0_linearModes( self ):

    lagrangeBasis = basisfunctions.Lagrange( )

    # r = -1
    linarModes = lagrangeBasis( -1.0, 1, 0 );
    self.assertEqual( linarModes.shape, ( 1, 2 ) )

    self.assertAlmostEqual( linarModes[0, 0], 1.0, delta = self.tolerance )
    self.assertAlmostEqual( linarModes[0, 1], 0.0, delta = self.tolerance )

    # r = 0
    linarModes = lagrangeBasis( 0.0, 1, 0 );
    self.assertEqual( linarModes.shape, ( 1, 2 ) )

    self.assertAlmostEqual( linarModes[0, 0], 0.5, delta = self.tolerance )
    self.assertAlmostEqual( linarModes[0, 1], 0.5, delta = self.tolerance )

    # r = 1
    linarModes = lagrangeBasis( 1.0, 1, 0 );
    self.assertEqual( linarModes.shape, ( 1, 2 ) )

    self.assertAlmostEqual( linarModes[0, 0], 0.0, delta = self.tolerance )
    self.assertAlmostEqual( linarModes[0, 1], 1.0, delta = self.tolerance )


  def test_difforder0_quadraticModes1( self ):

    lagrangeBasis = basisfunctions.Lagrange( )

    # r = -1
    quadraticModes = lagrangeBasis( -1.0, 2, 0 );
    self.assertEqual( quadraticModes.shape, ( 1, 3 ) )

    self.assertAlmostEqual( quadraticModes[0, 0], 1.0, delta = self.tolerance )
    self.assertAlmostEqual( quadraticModes[0, 1], 0.0, delta = self.tolerance )
    self.assertAlmostEqual( quadraticModes[0, 2], 0.0, delta = self.tolerance )

    # r = -0.5
    quadraticModes = lagrangeBasis( -0.5, 2, 0 );
    self.assertEqual( quadraticModes.shape, ( 1, 3 ) )

    self.assertAlmostEqual( quadraticModes[0, 0],  3.0 / 8.0, delta = self.tolerance )
    self.assertAlmostEqual( quadraticModes[0, 1],  3.0 / 4.0, delta = self.tolerance )
    self.assertAlmostEqual( quadraticModes[0, 2], -1.0 / 8.0, delta = self.tolerance )
  
    # r = 0
    quadraticModes = lagrangeBasis( 0.0, 2, 0 );
    self.assertEqual( quadraticModes.shape, ( 1, 3 ) )

    self.assertAlmostEqual( quadraticModes[0, 0], 0.0, delta = self.tolerance )
    self.assertAlmostEqual( quadraticModes[0, 1], 1.0, delta = self.tolerance )
    self.assertAlmostEqual( quadraticModes[0, 2], 0.0, delta = self.tolerance )
      
    # r = 0.5
    quadraticModes = lagrangeBasis( 0.5, 2, 0 );
    self.assertEqual( quadraticModes.shape, ( 1, 3 ) )

    self.assertAlmostEqual( quadraticModes[0, 0], -1.0 / 8.0, delta = self.tolerance )
    self.assertAlmostEqual( quadraticModes[0, 1],  3.0 / 4.0, delta = self.tolerance )
    self.assertAlmostEqual( quadraticModes[0, 2],  3.0 / 8.0, delta = self.tolerance )
  
    # r = 1
    quadraticModes = lagrangeBasis( 1.0, 2, 0 );
    self.assertEqual( quadraticModes.shape, ( 1, 3 ) )

    self.assertAlmostEqual( quadraticModes[0, 0], 0.0, delta = self.tolerance )
    self.assertAlmostEqual( quadraticModes[0, 1], 0.0, delta = self.tolerance )
    self.assertAlmostEqual( quadraticModes[0, 2], 1.0, delta = self.tolerance )

  def test_difforder0_quadraticModes2( self ):

    lagrangeBasis = basisfunctions.Lagrange( )

    analyticFunctions = [ lambda r: ( r * ( r - 1.0 ) ) / 2.0, 
                          lambda r: -( r - 1.0 ) * ( r + 1.0 ) , 
                          lambda r: r * ( r / 2.0 + 1.0 / 2.0 ) ];

    # choose where to evaluate
    coordinates = numpy.linspace( -1.0, 1.0, 13 );

    for i in xrange( coordinates.shape[0] ):

      shapes = lagrangeBasis( coordinates[i], 2, 0 );
      self.assertEqual( shapes.shape, ( 1, 3 ) )

      for j in xrange( len( analyticFunctions ) ):
        self.assertAlmostEqual( shapes[0, j], analyticFunctions[j]( coordinates[i] ), delta = self.tolerance )

  def test_difforder0_cubicModes( self ):

    lagrangeBasis = basisfunctions.Lagrange( )

    analyticFunctions = [ lambda r: -( 3.0 * ( ( 3.0 * r ) / 2.0 + 1.0 / 2.0 ) * ( r - 1.0 ) * ( r -  1.0 / 3.0 ) ) / 8.0, 
                          lambda r: ( 9.0 * ( ( 3.0 * r ) / 2.0 + 3.0 / 2.0 ) * ( r - 1.0 ) * ( r - 1.0 / 3.0 ) ) / 8.0 , 
                          lambda r: -( 9.0 * ( ( 3.0 * r ) / 4.0 + 3.0 / 4.0 ) * ( r - 1.0 ) * ( r + 1.0 / 3.0 ) ) / 4.0,
                          lambda r: ( 9.0 * ( r / 2.0 + 1.0 / 2.0 ) * ( r - 1.0 / 3.0 ) * ( r + 1.0 / 3.0 ) ) / 8.0 ];

    # choose where to evaluate
    coordinates = numpy.linspace( -1.0, 1.0, 13 );

    for i in xrange( coordinates.shape[0] ):

      shapes = lagrangeBasis( coordinates[i], 3, 0 );
      self.assertEqual( shapes.shape, ( 1, 4 ) )

      for j in xrange( len( analyticFunctions ) ):
        self.assertAlmostEqual( shapes[0, j], analyticFunctions[j]( coordinates[i] ), delta = self.tolerance )

  def test_difforder0_quarticModes( self ):

    lagrangeBasis = basisfunctions.Lagrange( )

    analyticFunctions = [ lambda r: ( r * ( 2.0 * r + 1.0 ) * ( r - 1.0 ) * ( r - 1.0 / 2.0 ) ) / 3.0, 
                          lambda r: -( 4.0 * r * ( 2.0 * r + 2.0 ) * ( r - 1.0 ) * ( r - 1.0 / 2.0 ) ) / 3.0, 
                          lambda r: 4.0 * ( r - 1.0 ) * ( r + 1.0 ) * ( r - 1.0 / 2.0 ) * ( r + 1.0 / 2.0 ),
                          lambda r: -4.0 * r * ( ( 2.0 * r ) / 3.0 + 2.0 / 3.0 ) * ( r - 1.0 ) * ( r + 1.0 / 2.0 ),
                          lambda r: ( 4.0 * r * ( r / 2.0  + 1.0 / 2.0 ) * ( r - 1.0 / 2.0 ) * ( r + 1.0 / 2.0 ) ) / 3.0 ];

    # choose where to evaluate
    coordinates = numpy.linspace( -1.0, 1.0, 13 );

    for i in xrange( coordinates.shape[0] ):

      shapes = lagrangeBasis( coordinates[i], 4, 0 );
      self.assertEqual( shapes.shape, ( 1, 5 ) )

      for j in xrange( len( analyticFunctions ) ):
        self.assertAlmostEqual( shapes[0, j], analyticFunctions[j]( coordinates[i] ), delta = self.tolerance )

  def test_difforder0_partitionOfUnity( self ):

    lagrangeBasis = basisfunctions.Lagrange( )

    polynomialDegrees = [ 3, 4, 5, 6, 8, 16 ]
    evaluationPoints = [ -0.66, -0.12, 0.23, 0.52, 0.76 ]

    # loop over all polynomial degrees and evaluation coordinates
    for p in polynomialDegrees:
      for r in evaluationPoints:

        # evaluate and check dimensions
        shapeFunctions = lagrangeBasis( r, p, 0 );
        self.assertEqual( shapeFunctions.shape, ( 1, p + 1 ) )

        # sum shape functions values up and check if sum is 1.0
        sumOfBasisFunctions = 0.0
        for i in xrange( 0, p + 1 ):
          sumOfBasisFunctions += shapeFunctions[0, i]

        self.assertAlmostEqual( sumOfBasisFunctions, 1.0, delta = self.tolerance )

      # for current p, check values at sample points i
      for i in xrange( 0, p + 1 ):          

        r = -1.0 + 2.0 * float(i) / float(p)

        # evaluate and check dimensions
        shapeFunctions = lagrangeBasis( r, p, 0 );
        self.assertEqual( shapeFunctions.shape, ( 1, p + 1 ) )

        # and loop over all basis functions
        for j in xrange( 0, p + 1 ):

          # if j == 1 shapes should be 1, else they should be zero
          self.assertAlmostEqual( shapeFunctions[0, j], float(i == j), delta = self.tolerance )

  def test_raiseException( self ):
    
    lagrangeBasis = basisfunctions.Lagrange( );

    self.assertRaises( ValueError, lagrangeBasis,  0.0,  2, -1 )

    # diff order 0
    self.assertRaises( ValueError, lagrangeBasis,  0.0, -1,  0 )
    self.assertRaises( ValueError, lagrangeBasis,  0.0,  0,  0 )
    self.assertRaises( ValueError, lagrangeBasis, -1.1,  1,  0 )
    self.assertRaises( ValueError, lagrangeBasis,  1.1,  0,  0 )

    # diff order 1
    self.assertRaises( ValueError, lagrangeBasis,  0.0, -1,  1 )
    self.assertRaises( ValueError, lagrangeBasis,  0.0,  0,  1 )
    self.assertRaises( ValueError, lagrangeBasis, -1.1,  1,  1 )
    self.assertRaises( ValueError, lagrangeBasis,  1.1,  0,  1 )

    # diff order 2
    self.assertRaises( ValueError, lagrangeBasis,  0.0, -1,  2 )
    self.assertRaises( ValueError, lagrangeBasis,  0.0,  0,  2 )
    self.assertRaises( ValueError, lagrangeBasis, -1.1,  1,  2 )
    self.assertRaises( ValueError, lagrangeBasis,  1.1,  0,  2 )

  def test_difforder1_linearModes( self ):

    lagrangeBasis = basisfunctions.Lagrange( )

    # r = -1
    linarModes = lagrangeBasis( -0.7, 1, 1 );
    self.assertEqual( linarModes.shape, ( 1, 2 ) )

    self.assertAlmostEqual( linarModes[0, 0], -0.5, delta = self.tolerance )
    self.assertAlmostEqual( linarModes[0, 1],  0.5, delta = self.tolerance )

  def test_difforder1_quadraticModes( self ):

    lagrangeBasis = basisfunctions.Lagrange( )

    analyticFunctions = [ lambda r: r - 1.0 / 2.0, 
                          lambda r: -2.0 * r, 
                          lambda r: r + 1.0 / 2.0 ];

    # choose where to evaluate
    coordinates = numpy.linspace( -1.0, 1.0, 13 );

    for i in xrange( coordinates.shape[0] ):

      shapes = lagrangeBasis( coordinates[i], 2, 1 );
      self.assertEqual( shapes.shape, ( 1, 3 ) )

      for j in xrange( len( analyticFunctions ) ):
        self.assertAlmostEqual( shapes[0, j], analyticFunctions[j]( coordinates[i] ), delta = self.tolerance )

  def test_difforder1_cubicModes( self ):

    lagrangeBasis = basisfunctions.Lagrange( )

    analyticFunctions = [ lambda r: ( 9.0 * r ) / 8.0 - ( 27.0 * r * r ) / 16.0 + 1.0 / 16.0, 
                          lambda r: ( 81.0* r * r ) / 16.0 - ( 9.0 * r ) / 8.0 - 27.0 / 16.0, 
                          lambda r: 27.0 / 16.0 - ( 81.0 * r * r ) / 16.0 - ( 9.0 * r ) / 8.0,
                          lambda r: ( 27.0 * r * r ) / 16.0 + ( 9.0 * r ) / 8.0 - 1.0 / 16.0 ];

    # choose where to evaluate
    coordinates = numpy.linspace( -1.0, 1.0, 13 );

    for i in xrange( coordinates.shape[0] ):

      shapes = lagrangeBasis( coordinates[i], 3, 1 );
      self.assertEqual( shapes.shape, ( 1, 4 ) )

      for j in xrange( len( analyticFunctions ) ):
        self.assertAlmostEqual( shapes[0, j], analyticFunctions[j]( coordinates[i] ), delta = self.tolerance )

  def test_difforder1_quarticModes( self ):

    lagrangeBasis = basisfunctions.Lagrange( )

    analyticFunctions = [ lambda r: ( 8.0 * r * r * r ) / 3.0 - 2.0 * r * r - r / 3.0 + 1.0 / 6.0, 
                          lambda r: -( 32.0 * r * r * r ) / 3.0 + 4.0 * r * r  + ( 16.0 * r ) / 3.0 - 4.0 / 3.0, 
                          lambda r: 2.0 * r * ( 8.0 * r * r - 5.0 ),
                          lambda r: -( 32.0 * r * r * r ) / 3.0 - 4.0 * r * r + ( 16.0 * r ) / 3.0 + 4.0 / 3.0,
                          lambda r: ( 8.0 * r * r * r ) / 3.0 + 2.0 * r * r - r / 3.0 - 1.0 / 6.0 ];

    # choose where to evaluate
    coordinates = numpy.linspace( -1.0, 1.0, 13 );

    for i in xrange( coordinates.shape[0] ):

      shapes = lagrangeBasis( coordinates[i], 4, 1 );
      self.assertEqual( shapes.shape, ( 1, 5 ) )

      for j in xrange( len( analyticFunctions ) ):
        self.assertAlmostEqual( shapes[0, j], analyticFunctions[j]( coordinates[i] ), delta = self.tolerance )

  def test_difforder2_linearcModes( self ):

    lagrangeBasis = basisfunctions.Lagrange( )

    analyticFunctions = [ lambda r: 0.0, 
                          lambda r: 0.0 ];

    # choose where to evaluate
    coordinates = numpy.linspace( -1.0, 1.0, 13 );

    for i in xrange( coordinates.shape[0] ):

      shapes = lagrangeBasis( coordinates[i], 1, 2 );
      self.assertEqual( shapes.shape, ( 1, 2 ) )

      for j in xrange( len( analyticFunctions ) ):
        self.assertAlmostEqual( shapes[0, j], analyticFunctions[j]( coordinates[i] ), delta = self.tolerance )

  def test_difforder2_quadraticModes( self ):

    lagrangeBasis = basisfunctions.Lagrange( )

    analyticFunctions = [ lambda r: 1.0, 
                          lambda r: -2.0, 
                          lambda r: 1.0 ];

    # choose where to evaluate
    coordinates = numpy.linspace( -1.0, 1.0, 13 );

    for i in xrange( coordinates.shape[0] ):

      shapes = lagrangeBasis( coordinates[i], 2, 2 );
      self.assertEqual( shapes.shape, ( 1, 3 ) )

      for j in xrange( len( analyticFunctions ) ):
        self.assertAlmostEqual( shapes[0, j], analyticFunctions[j]( coordinates[i] ), delta = self.tolerance )

  def test_difforder2_cubicModes( self ):

    lagrangeBasis = basisfunctions.Lagrange( )

    analyticFunctions = [ lambda r: 9.0 / 8.0 - ( 27.0 * r ) / 8.0, 
                          lambda r: (81.0 * r ) / 8.0 - 9.0 / 8.0, 
                          lambda r: -( 81.0 * r ) / 8.0 - 9.0 / 8.0, 
                          lambda r: ( 27.0 * r ) / 8.0 + 9.0 / 8.0 ];

    # choose where to evaluate
    coordinates = numpy.linspace( -1.0, 1.0, 13 );

    for i in xrange( coordinates.shape[0] ):

      shapes = lagrangeBasis( coordinates[i], 3, 2 );
      self.assertEqual( shapes.shape, ( 1, 4 ) )

      for j in xrange( len( analyticFunctions ) ):
        self.assertAlmostEqual( shapes[0, j], analyticFunctions[j]( coordinates[i] ), delta = self.tolerance )

  def test_difforder2_quarticModes( self ):

    lagrangeBasis = basisfunctions.Lagrange( )

    analyticFunctions = [ lambda r: 8.0 * r * r - 4.0 * r - 1.0 / 3.0, 
                          lambda r: -32.0 * r * r + 8.0 * r + 16.0 / 3.0, 
                          lambda r: 48.0 * r * r - 10.0, 
                          lambda r: -32.0 * r * r - 8.0 * r + 16.0 / 3.0, 
                          lambda r: 8.0 * r * r + 4.0 * r - 1.0 / 3.0 ];

    # choose where to evaluate
    coordinates = numpy.linspace( -1.0, 1.0, 13 );

    for i in xrange( coordinates.shape[0] ):

      shapes = lagrangeBasis( coordinates[i], 4, 2 );
      self.assertEqual( shapes.shape, ( 1, 5 ) )

      for j in xrange( len( analyticFunctions ) ):
        self.assertAlmostEqual( shapes[0, j], analyticFunctions[j]( coordinates[i] ), delta = self.tolerance )


class TestGridTensorProduct(unittest.TestCase):

  def setUp( self ):
    self.tolerance = 1e-10

  def test_LinearShapes2D( self ):

    for basis1D in [ basisfunctions.IntegratedLegendre( ), basisfunctions.Lagrange( ) ]:

      tensorProduct = basisfunctions.GridTensorProduct( basis1D ) 

      rst = [ [ -1.0, -0.43, 0.0, 0.14, 1.0 ], [ -1.0, -0.87, -0.04, 0.0, 0.61, 0.89, 1.0 ] ];

      shapes = tensorProduct.evaluate( rst, 1, 2, 3 )

      for iEvaluationPoint in xrange( len(rst[0]) ):
        for jEvaluationPoint in xrange( len(rst[1]) ):
          r = rst[0][iEvaluationPoint]
          s = rst[1][jEvaluationPoint]

          for iFieldComponent in xrange( 3 ):

            N = shapes[iFieldComponent][0][iEvaluationPoint, jEvaluationPoint];

            dNdx = shapes[iFieldComponent][1][iEvaluationPoint, jEvaluationPoint, 0];
            dNdy = shapes[iFieldComponent][1][iEvaluationPoint, jEvaluationPoint, 1];

            ddNdxdx = shapes[iFieldComponent][2][iEvaluationPoint, jEvaluationPoint, 0];
            ddNdydy = shapes[iFieldComponent][2][iEvaluationPoint, jEvaluationPoint, 1];
            ddNdxdy = shapes[iFieldComponent][2][iEvaluationPoint, jEvaluationPoint, 2];

            # diff order 0
            self.assertAlmostEqual( N[0], 0.25 * ( 1.0 - r ) * ( 1.0 - s ), delta = self.tolerance )
            self.assertAlmostEqual( N[1], 0.25 * ( 1.0 - r ) * ( 1.0 + s ), delta = self.tolerance )
            self.assertAlmostEqual( N[2], 0.25 * ( 1.0 + r ) * ( 1.0 - s ), delta = self.tolerance )
            self.assertAlmostEqual( N[3], 0.25 * ( 1.0 + r ) * ( 1.0 + s ), delta = self.tolerance )
                       
            # diff order 1
            self.assertAlmostEqual( dNdx[0], 0.25 * ( - 1.0 ) * ( 1.0 - s ), delta = self.tolerance )
            self.assertAlmostEqual( dNdx[1], 0.25 * ( - 1.0 ) * ( 1.0 + s ), delta = self.tolerance )
            self.assertAlmostEqual( dNdx[2], 0.25 * (   1.0 ) * ( 1.0 - s ), delta = self.tolerance )
            self.assertAlmostEqual( dNdx[3], 0.25 * (   1.0 ) * ( 1.0 + s ), delta = self.tolerance )
                                   
            self.assertAlmostEqual( dNdy[0], 0.25 * ( 1.0 - r ) * ( -1.0 ), delta = self.tolerance )
            self.assertAlmostEqual( dNdy[1], 0.25 * ( 1.0 - r ) * (  1.0 ), delta = self.tolerance )
            self.assertAlmostEqual( dNdy[2], 0.25 * ( 1.0 + r ) * ( -1.0 ), delta = self.tolerance )
            self.assertAlmostEqual( dNdy[3], 0.25 * ( 1.0 + r ) * (  1.0 ), delta = self.tolerance )

            # diff order 2:
            self.assertAlmostEqual( ddNdxdx[0], 0.0, delta = self.tolerance )
            self.assertAlmostEqual( ddNdxdx[1], 0.0, delta = self.tolerance )
            self.assertAlmostEqual( ddNdxdx[2], 0.0, delta = self.tolerance )
            self.assertAlmostEqual( ddNdxdx[3], 0.0, delta = self.tolerance )
                                   
            self.assertAlmostEqual( ddNdydy[0], 0.0, delta = self.tolerance )
            self.assertAlmostEqual( ddNdydy[1], 0.0, delta = self.tolerance )
            self.assertAlmostEqual( ddNdydy[2], 0.0, delta = self.tolerance )
            self.assertAlmostEqual( ddNdydy[3], 0.0, delta = self.tolerance )
                                                     
            self.assertAlmostEqual( ddNdxdy[0], 0.25 * ( -1.0 ) * ( -1.0 ), delta = self.tolerance )
            self.assertAlmostEqual( ddNdxdy[1], 0.25 * ( -1.0 ) * (  1.0 ), delta = self.tolerance )
            self.assertAlmostEqual( ddNdxdy[2], 0.25 * (  1.0 ) * ( -1.0 ), delta = self.tolerance )
            self.assertAlmostEqual( ddNdxdy[3], 0.25 * (  1.0 ) * (  1.0 ), delta = self.tolerance )
                  

  def test_LinearShapes2DAnisotropLagrange( self ):

    basis1D = basisfunctions.Lagrange( )
    tensorProduct = basisfunctions.GridTensorProduct( basis1D ) 

    # 7 points in r, 13 points in s
    rst = [ numpy.linspace(-1.0, 1.0, 7), 
            numpy.linspace(-1.0, 1.0, 13) ];

    #rst = [ numpy.linspace(-1.0, 1.0, 3), 
    #        [-1.0] ];

    # 3 field components, 2 space dimensions
    p = [ [ 3, 2 ],
          [ 2, 3 ], 
          [ 2, 2 ] ];

    shapes = tensorProduct.evaluate( rst, p, 2, 3 )

    # Quadratic Lagrange polynomials up to second derivative
    quadraticLagrange = [ [ lambda r: ( r * ( r - 1.0 ) ) / 2.0, 
                            lambda r: -( r - 1.0 ) * ( r + 1.0 ) , 
                            lambda r: r * ( r / 2.0 + 1.0 / 2.0 ) ],
                          [ lambda r: r - 1.0 / 2.0, 
                            lambda r: -2.0 * r, 
                            lambda r: r + 1.0 / 2.0 ],
                          [ lambda r: 1.0, 
                            lambda r: -2.0, 
                            lambda r: 1.0 ] ]

    # Cubic Lagrange polynomials up to second derivative
    cubicLagrange = [ [ lambda r: -( 3.0 * ( ( 3.0 * r ) / 2.0 + 1.0 / 2.0 ) * ( r - 1.0 ) * ( r -  1.0 / 3.0 ) ) / 8.0, 
                        lambda r: ( 9.0 * ( ( 3.0 * r ) / 2.0 + 3.0 / 2.0 ) * ( r - 1.0 ) * ( r - 1.0 / 3.0 ) ) / 8.0 , 
                        lambda r: -( 9.0 * ( ( 3.0 * r ) / 4.0 + 3.0 / 4.0 ) * ( r - 1.0 ) * ( r + 1.0 / 3.0 ) ) / 4.0,
                        lambda r: ( 9.0 * ( r / 2.0 + 1.0 / 2.0 ) * ( r - 1.0 / 3.0 ) * ( r + 1.0 / 3.0 ) ) / 8.0 ],
                      [ lambda r: ( 9.0 * r ) / 8.0 - ( 27.0 * r * r ) / 16.0 + 1.0 / 16.0, 
                        lambda r: ( 81.0* r * r ) / 16.0 - ( 9.0 * r ) / 8.0 - 27.0 / 16.0, 
                        lambda r: 27.0 / 16.0 - ( 81.0 * r * r ) / 16.0 - ( 9.0 * r ) / 8.0,
                        lambda r: ( 27.0 * r * r ) / 16.0 + ( 9.0 * r ) / 8.0 - 1.0 / 16.0 ],
                      [ lambda r: 9.0 / 8.0 - ( 27.0 * r ) / 8.0, 
                        lambda r: (81.0 * r ) / 8.0 - 9.0 / 8.0, 
                        lambda r: -( 81.0 * r ) / 8.0 - 9.0 / 8.0, 
                        lambda r: ( 27.0 * r ) / 8.0 + 9.0 / 8.0 ] ]

    for iEvaluationPoint in xrange( len(rst[0]) ):
      for jEvaluationPoint in xrange( len(rst[1]) ):

        r = rst[0][iEvaluationPoint]
        s = rst[1][jEvaluationPoint]

        # ------------------------ field component 0 --------------------------          
        N = shapes[0][0][iEvaluationPoint, jEvaluationPoint];

        dNdx = shapes[0][1][iEvaluationPoint, jEvaluationPoint, 0];
        dNdy = shapes[0][1][iEvaluationPoint, jEvaluationPoint, 1];

        ddNdxdx = shapes[0][2][iEvaluationPoint, jEvaluationPoint, 0];
        ddNdydy = shapes[0][2][iEvaluationPoint, jEvaluationPoint, 1];
        ddNdxdy = shapes[0][2][iEvaluationPoint, jEvaluationPoint, 2];

        def testTensorProduct( diff1, diff2, shapes ):

          self.assertEqual( len(shapes), 12 )

          self.assertAlmostEqual( shapes[ 0], cubicLagrange[diff1][0](r) * quadraticLagrange[diff2][0](s), delta=self.tolerance )
          self.assertAlmostEqual( shapes[ 1], cubicLagrange[diff1][0](r) * quadraticLagrange[diff2][1](s), delta=self.tolerance )
          self.assertAlmostEqual( shapes[ 2], cubicLagrange[diff1][0](r) * quadraticLagrange[diff2][2](s), delta=self.tolerance )
          self.assertAlmostEqual( shapes[ 3], cubicLagrange[diff1][1](r) * quadraticLagrange[diff2][0](s), delta=self.tolerance )
          self.assertAlmostEqual( shapes[ 4], cubicLagrange[diff1][1](r) * quadraticLagrange[diff2][1](s), delta=self.tolerance )
          self.assertAlmostEqual( shapes[ 5], cubicLagrange[diff1][1](r) * quadraticLagrange[diff2][2](s), delta=self.tolerance )
          self.assertAlmostEqual( shapes[ 6], cubicLagrange[diff1][2](r) * quadraticLagrange[diff2][0](s), delta=self.tolerance )
          self.assertAlmostEqual( shapes[ 7], cubicLagrange[diff1][2](r) * quadraticLagrange[diff2][1](s), delta=self.tolerance )
          self.assertAlmostEqual( shapes[ 8], cubicLagrange[diff1][2](r) * quadraticLagrange[diff2][2](s), delta=self.tolerance )
          self.assertAlmostEqual( shapes[ 9], cubicLagrange[diff1][3](r) * quadraticLagrange[diff2][0](s), delta=self.tolerance )
          self.assertAlmostEqual( shapes[10], cubicLagrange[diff1][3](r) * quadraticLagrange[diff2][1](s), delta=self.tolerance )
          self.assertAlmostEqual( shapes[11], cubicLagrange[diff1][3](r) * quadraticLagrange[diff2][2](s), delta=self.tolerance )

        testTensorProduct( 0, 0, N )
        testTensorProduct( 1, 0, dNdx )
        testTensorProduct( 0, 1, dNdy )
        testTensorProduct( 2, 0, ddNdxdx )
        testTensorProduct( 0, 2, ddNdydy )
        testTensorProduct( 1, 1, ddNdxdy )

        # ------------------------ field component 1 --------------------------          
        N = shapes[1][0][iEvaluationPoint, jEvaluationPoint];

        dNdx = shapes[1][1][iEvaluationPoint, jEvaluationPoint, 0];
        dNdy = shapes[1][1][iEvaluationPoint, jEvaluationPoint, 1];

        ddNdxdx = shapes[1][2][iEvaluationPoint, jEvaluationPoint, 0];
        ddNdydy = shapes[1][2][iEvaluationPoint, jEvaluationPoint, 1];
        ddNdxdy = shapes[1][2][iEvaluationPoint, jEvaluationPoint, 2];

        def testTensorProduct( diff1, diff2, shapes ):

          self.assertEqual( len(shapes), 12 )

          self.assertAlmostEqual( shapes[ 0], quadraticLagrange[diff1][0](r) * cubicLagrange[diff2][0](s), delta=self.tolerance )
          self.assertAlmostEqual( shapes[ 1], quadraticLagrange[diff1][0](r) * cubicLagrange[diff2][1](s), delta=self.tolerance )
          self.assertAlmostEqual( shapes[ 2], quadraticLagrange[diff1][0](r) * cubicLagrange[diff2][2](s), delta=self.tolerance )
          self.assertAlmostEqual( shapes[ 3], quadraticLagrange[diff1][0](r) * cubicLagrange[diff2][3](s), delta=self.tolerance )
          self.assertAlmostEqual( shapes[ 4], quadraticLagrange[diff1][1](r) * cubicLagrange[diff2][0](s), delta=self.tolerance )
          self.assertAlmostEqual( shapes[ 5], quadraticLagrange[diff1][1](r) * cubicLagrange[diff2][1](s), delta=self.tolerance )
          self.assertAlmostEqual( shapes[ 6], quadraticLagrange[diff1][1](r) * cubicLagrange[diff2][2](s), delta=self.tolerance )
          self.assertAlmostEqual( shapes[ 7], quadraticLagrange[diff1][1](r) * cubicLagrange[diff2][3](s), delta=self.tolerance )
          self.assertAlmostEqual( shapes[ 8], quadraticLagrange[diff1][2](r) * cubicLagrange[diff2][0](s), delta=self.tolerance )
          self.assertAlmostEqual( shapes[ 9], quadraticLagrange[diff1][2](r) * cubicLagrange[diff2][1](s), delta=self.tolerance )
          self.assertAlmostEqual( shapes[10], quadraticLagrange[diff1][2](r) * cubicLagrange[diff2][2](s), delta=self.tolerance )
          self.assertAlmostEqual( shapes[11], quadraticLagrange[diff1][2](r) * cubicLagrange[diff2][3](s), delta=self.tolerance )

        testTensorProduct( 0, 0, N )
        testTensorProduct( 1, 0, dNdx )
        testTensorProduct( 0, 1, dNdy )
        testTensorProduct( 2, 0, ddNdxdx )
        testTensorProduct( 0, 2, ddNdydy )
        testTensorProduct( 1, 1, ddNdxdy )

        # ------------------------ field component 2 --------------------------          
        N = shapes[2][0][iEvaluationPoint, jEvaluationPoint];

        dNdx = shapes[2][1][iEvaluationPoint, jEvaluationPoint, 0];
        dNdy = shapes[2][1][iEvaluationPoint, jEvaluationPoint, 1];

        ddNdxdx = shapes[2][2][iEvaluationPoint, jEvaluationPoint, 0];
        ddNdydy = shapes[2][2][iEvaluationPoint, jEvaluationPoint, 1];
        ddNdxdy = shapes[2][2][iEvaluationPoint, jEvaluationPoint, 2];

        def testTensorProduct( diff1, diff2, shapes ):

          self.assertEqual( len(shapes), 9 )

          self.assertAlmostEqual( shapes[0], quadraticLagrange[diff1][0](r) * quadraticLagrange[diff2][0](s), delta=self.tolerance )
          self.assertAlmostEqual( shapes[1], quadraticLagrange[diff1][0](r) * quadraticLagrange[diff2][1](s), delta=self.tolerance )
          self.assertAlmostEqual( shapes[2], quadraticLagrange[diff1][0](r) * quadraticLagrange[diff2][2](s), delta=self.tolerance )
          self.assertAlmostEqual( shapes[3], quadraticLagrange[diff1][1](r) * quadraticLagrange[diff2][0](s), delta=self.tolerance )
          self.assertAlmostEqual( shapes[4], quadraticLagrange[diff1][1](r) * quadraticLagrange[diff2][1](s), delta=self.tolerance )
          self.assertAlmostEqual( shapes[5], quadraticLagrange[diff1][1](r) * quadraticLagrange[diff2][2](s), delta=self.tolerance )
          self.assertAlmostEqual( shapes[6], quadraticLagrange[diff1][2](r) * quadraticLagrange[diff2][0](s), delta=self.tolerance )
          self.assertAlmostEqual( shapes[7], quadraticLagrange[diff1][2](r) * quadraticLagrange[diff2][1](s), delta=self.tolerance )
          self.assertAlmostEqual( shapes[8], quadraticLagrange[diff1][2](r) * quadraticLagrange[diff2][2](s), delta=self.tolerance )

        testTensorProduct( 0, 0, N )
        testTensorProduct( 1, 0, dNdx )
        testTensorProduct( 0, 1, dNdy )
        testTensorProduct( 2, 0, ddNdxdx )
        testTensorProduct( 0, 2, ddNdydy )
        testTensorProduct( 1, 1, ddNdxdy )

if __name__ == '__main__':
    unittest.main()







