# --- Python Imports ---
import numpy as np

# ---------------------------------------------------------

class Integrator:
    '''
    Gaussian quadrature
    '''
    def __init__( self, polynomialOrder ):
        if polynomialOrder < 1:
            raise AttributeError( "Invalid polynomial order!" )

        # Get sample points and weights
        self.polynomialOrder = polynomialOrder
        self.positions, self.weights = np.polynomial.legendre.leggauss( polynomialOrder )



    def __call__( self, scalarFunctor, domain ):
        # Check domain
        scale = (float(domain[1]) - float(domain[0])) / 2.0
        if scale <= 0:
            raise ValueError( "Invalid integration domain!" )

        # Integrate
        value = 0.0
        for position, weight in zip( self.positions, self.weights ):
            value += weight * scalarFunctor( scale*(position+1.0) + domain[0] )
        
        return value*scale


    def createCache( self, scalarFunctor, domain ):
        # Check domain
        scale = (float(domain[1]) - float(domain[0])) / 2.0
        if scale <= 0:
            raise ValueError( "Invalid integration domain!" )

        return np.asarray([ scalarFunctor( scale*(position+1.0) + domain[0] ) for position in self.positions ])


    def integrateCached( self, scalarFunctor, cachedCoefficients, domain ):
        # Check domain
        scale = (float(domain[1]) - float(domain[0])) / 2.0
        if scale <= 0:
            raise ValueError( "Invalid integration domain!" )

        # Integrate
        value = 0.0
        for position, weight, cache in zip( self.positions, self.weights, cachedCoefficients ):
            value += weight * cache * scalarFunctor( scale*(position+1.0) + domain[0] )
        
        return value*scale