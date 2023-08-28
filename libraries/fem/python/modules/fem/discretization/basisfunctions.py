# --- Python Imports ---
import numpy as np
from numpy.polynomial.legendre import Legendre, legint

# --- Utility Imports ---
from cie.fem.utilities import *

# ---------------------------------------------------------
class BasisFunctions:
    '''
    Abstract base class for a matching collection of 1D basis functions

    Constructor arguments:
        functions=None          : basis functions to store
        dtype=np.float64        : storage and function return data types
    '''
    def __init__( self, functions=None, dtype=np.float64 ):
        # Initialize basic members
        self.dtype              = dtype
        self.domain             = np.array( (-1.0, 1.0), dtype=dtype )
        self._functions         = functions
        self._nodalFunctions    = None
        self._counter           = -1


    def evaluate( self, basisID, position ):
        # Check if the function was set
        if self._functions == None:
            raise AttributeError("Unset basis function!")

        # Define filter
        def windowedBasis( pos ):
            if self.domain[0]<=pos and pos<=self.domain[1]:
                return self.dtype( self._functions[basisID](pos) )
            else:
                return self.dtype(0.0)

        # Evaluate at a single point
        if not isNumpyArray(position):
            return windowedBasis(position)

        # Evaluate at multiple points
        else:
            return np.asarray( [ windowedBasis(x) for x in position], dtype=self.dtype )



    def evaluateAll( self, position ):
        if isArray( position ):
            return np.asarray( [self.evaluate(basisID,position) for basisID in range(len(self._functions))] )
        else:
            raise ValueError( "Expecting an array argument!" )


    def __call__( self, basisID, position ):
        return self.evaluate( basisID, position )


    def __iter__( self ):
        self._counter = -1
        return self


    def __next__( self ):
        if self._counter < len(self._functions)-1:
            self._counter += 1
            return self._functions[self._counter]
        else:
            raise StopIteration


    def __len__( self ):
        return len(self._functions)


    def __getitem__( self, index ):
        return self._functions[index]


    @property
    def nodalFunctions(self):
        # Check whether finding nodal functions == necessary
        if self._nodalFunctions != None:
            return self._nodalFunctions

        else:
            # Check whether instance was initialized
            if self._functions == None:
                raise AttributeError( "List of basis functions must be filled before calling initialize" )

            # Find functions that do not vanish at the borders of the domain
            # (needed during assembly, to ensure C0 continuity across elements)
            self._nodalFunctions   = []
            for i, function in enumerate(self._functions):
                if np.abs(function(self.domain[0]))>1e-12 or np.abs(function(self.domain[1]))>1e-12:
                    self._nodalFunctions.append(i)

            return self._nodalFunctions


class PolynomialBasisFunctions( BasisFunctions ):
    '''
    1D polynomial basis function abstract base class

    Constructor arguments:
        polynomialOrder=1       : polynomial order of the basis to compute the number of functions
        functions=None          : basis functions to store
        dtype=np.float64        : storage and function return data types
    '''
    def __init__( self, *args, polynomialOrder=1, **kwargs ):
        BasisFunctions.__init__( self, *args, **kwargs )
        self.polynomialOrder = polynomialOrder



class LinearBasisFunctions( PolynomialBasisFunctions ):
    '''
    Classic linear basis functions

    Constructor arguments:
        dtype=np.float64        : storage and function return data types
    '''
    def __init__( self, *args, **kwargs ):
        PolynomialBasisFunctions.__init__( self, *args, **kwargs )
        self._functions = [ lambda t: 0.5 - t/2.0,
                            lambda t: 0.5 + t/2.0 ]


    def derivatives( self ):
        return BasisFunctions(  functions = [   lambda t: -0.5,
                                                lambda t: 0.5 ] )



class HierarchicBasisFunctions( PolynomialBasisFunctions ):
    '''
    Legendre polynomials with the constant function replaced by a linear one
    (using the numpy implementation of Legendre polynomials)

    Constructor arguments:
        polynomialOrder=1       : polynomial order of the basis to compute the number of functions
        dtype=np.float64        : storage and function return data types
    '''
    def __init__( self, *args, **kwargs ):
        PolynomialBasisFunctions.__init__( self, *args, **kwargs )

        # Standard linear basis functions if the polynomial order == 1
        if self.polynomialOrder == 1:
            self._functions = [ Legendre( (0.5,-0.5) ),
                                Legendre( (0.5,0.5) ) ]

        # Legendre polynomials, but replace the constant function with a linear one
        elif self.polynomialOrder > 1:
            self._functions = [None for i in range(self.polynomialOrder+1)]
            for index, coefficients in enumerate( oneHotArray(self.polynomialOrder+1) ):

                if index == 0:
                    coefficients[0] = 0
                    coefficients[1] = -1

                self._functions[index] = Legendre(  coefficients,
                                                    domain=self.domain  )

        else:
            raise ValueError( "Invalid polynomial order!" )


    def derivatives( self ):
        functions = [ Legendre(np.polynomial.legendre.legder(polynomial.coef)) for polynomial in self._functions ]
        return BasisFunctions( functions=functions )




class IntegratedHierarchicBasisFunctions( PolynomialBasisFunctions ):
    '''
    Hierarchic basis functions as proposed by Szabó & Babuska
    For easier assembly, the first and last functions are nodal ones, the others vanish at the boundaries

    Constructor arguments:
        polynomialOrder=1       : polynomial order of the basis to compute the number of functions
        dtype=np.float64        : storage and function return data types
    '''
    def __init__( self, *args, **kwargs ):
        PolynomialBasisFunctions.__init__( self, *args, **kwargs )

        # Standard linear basis functions if the polynomial order == 1
        if self.polynomialOrder == 1:
            self._functions = [ Legendre( (0.5,-0.5) ),
                                Legendre( (0.5,0.5) ) ]

        # Integrated Legendre polynomials (first two functions replaced by classic linear ones)
        elif self.polynomialOrder > 1:
            # Base
            self._functions     = [ Legendre( coefficients ) for coefficients in oneHotArray(self.polynomialOrder+1) ]

            # Integrate base functions
            for index in range( self.polynomialOrder-1, 0, -1 ):
                self._functions[index] = Legendre( legint( np.sqrt(float(index) - 0.5) * self._functions[index].coef,
                                                                    m=1,
                                                                    lbnd=-1.0),
                                                            domain=self.domain  )

            # Replace first and last functions
            self._functions[0]  = Legendre( (0.5,-0.5) )
            self._functions[-1] = Legendre( (0.5, 0.5) )

        else:
            raise ValueError( "Invalid polynomial order!" )


    def derivatives( self ):
        functions = [ Legendre(np.polynomial.legendre.legder(polynomial.coef)) for polynomial in self._functions ]
        return BasisFunctions( functions=functions )

