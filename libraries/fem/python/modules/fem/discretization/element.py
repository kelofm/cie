# --- Python Imports ---
import numpy as np
import collections

# --- Internal Imports ---
from .basisfunctions import *
from cie.fem.numeric import Integrator
from cie.fem.utilities import isArray
from cie.fem.utilities.cache import Cache 

# ---------------------------------------------------------

class Element:
    '''
    Element base class
    '''
    def __init__( self, domain, DoFs, load, basisFunctions=LinearBasisFunctions(), integrationOrder=2 ):
        # Check degrees of freedom
        if len(DoFs) is not len(basisFunctions._functions):
            raise AttributeError( "Number of degrees of freedom must match" )

        # Set members
        self.basisFunctions         = basisFunctions
        self.domain                 = np.asarray( domain )
        self.DoFs                   = np.asarray( DoFs )
        self.load                   = load
        self.integrator             = Integrator(integrationOrder)

        self.basisDerivatives       = basisFunctions.derivatives()

        self._basisDerivativeCache  = None
        self._basisCache            = None
        self._sampleCache           = Cache()
        self._derivativeCache       = Cache()

    
    def copy( self ):
        copy = Element( self.domain,
                        self.DoFs,
                        self.load,
                        basisFunctions=self.basisFunctions,
                        integrationOrder=self.integrator.polynomialOrder )

        # Do not copy the integrator and its caches
        return copy


    def updateGlobalMatrix( self, matrix, index1, index2, value ):
        matrix[self.DoFs[index1], self.DoFs[index2]] += value
        if index1 != index2:
            matrix[self.DoFs[index2], self.DoFs[index1]] += value


    def updateGlobalVector( self, vector, index, value ):
        vector[self.DoFs[index]] += value


    def integrateStiffness( self, globalStiffnessMatrix ):
        raise NotImplementedError( "Attempt to call abstract base class function" )


    def integrateMass( self, globalMassMatrix ):
        raise NotImplementedError( "Attempt to call abstract base class function" )


    def integrateLoad( self, globalLoadVector ):
        raise NotImplementedError( "Attempt to call abstract base class function" )



class Element1D( Element ):
    '''
    '''
    def __init__( self, *args, **kwargs ):
        Element.__init__( self, *args, **kwargs )
        
        # Check domain
        if self.domain[1] <= self.domain[0]:
            self.domain = np.flip( self.domain )

        # Jacobian ( dx/dxi )
        self._jacobian      = (self.domain[1]-self.domain[0])/2.0
        self._invJacobian   = 1.0/self._jacobian

    
    def copy( self ):
        copy = Element1D(   self.domain,
                            self.DoFs,
                            self.load,
                            basisFunctions=self.basisFunctions,
                            integrationOrder=self.integrator.polynomialOrder)
        return copy


    def toGlobalCoordinates( self, localCoordinates ):
        '''
        [-1,1] -> self.domain
        '''
        if isArray( localCoordinates ):
            return np.asarray( [ self._jacobian * (xi + 1.0) + self.domain[0] for xi in localCoordinates ] )
        else:
            return self._jacobian * (localCoordinates + 1.0) + self.domain[0]


    def toLocalCoordinates( self, globalCoordinates ):
        '''
        self.domain -> [-1,1]
        '''
        if isArray( globalCoordinates ):
            return np.asarray( [ (x-self.domain[0]) * self._invJacobian - 1.0 for x in globalCoordinates ] )
        else:
            return (globalCoordinates-self.domain[0]) * self._invJacobian - 1.0


    def __call__( self, coefficients, positions ):
        '''
        Given the solution coefficients, sample and sum the basis functions at the given coordinates
        '''
        # Initialize
        values = None
        if isNumpyArray( positions ):
            values = np.zeros( positions.shape, positions.dtype )
        else:
            values = 0.0

        # Get basis functions
        cacheID = self._sampleCache.hash(positions)
        if not self._sampleCache.check(cacheID, hashed=True):
            localCoordinates = self.toLocalCoordinates(positions)
            self._sampleCache.overwrite( cacheID, [self.basisFunctions( basisID, localCoordinates ) for basisID in range(len(self.basisFunctions)) ], hashed=True )
        

        # Compute solution
        for coefficient, basisValues in zip( coefficients, self._sampleCache.get(cacheID, hashed=True) ):
            #values += coefficient * self.basisFunctions( basisID, self.toLocalCoordinates(positions) )
            values += coefficient * basisValues

        return values


    def derivative( self, coefficients, positions ):
        '''
        Given the solution coefficients, sample and sum the derivative 
        of the basis functions at the given coordinates
        '''
        # Initialize
        values = None
        if isNumpyArray( positions ):
            values = np.zeros( positions.shape, positions.dtype )
        else:
            values = 0.0

        # Get basis functions
        cacheID = self._derivativeCache.hash(positions)
        if not self._derivativeCache.check(cacheID, hashed=True):
            localCoordinates = self.toLocalCoordinates(positions)
            self._derivativeCache.overwrite( cacheID, [self._invJacobian * self.basisDerivatives( basisID, localCoordinates ) for basisID in range(len(self.basisDerivatives)) ], hashed=True )
        

        # Compute solution
        for coefficient, derivativeValues in zip( coefficients, self._derivativeCache.get(cacheID, hashed=True) ):
            values += coefficient * derivativeValues

        return values


class LinearHeatElement1D( Element1D ):
    '''
    Most basic element for the Poisson equation
    '''
    def __init__( self, capacity, conductivity, *args, **kwargs ):
        Element1D.__init__(self, *args, **kwargs )
        self.capacity       = capacity
        self.conductivity   = conductivity


    def copy( self ):
        copy = LinearHeatElement1D( self.capacity,
                                    self.conductivity,
                                    self.domain,
                                    self.DoFs,
                                    self.load,
                                    basisFunctions=self.basisFunctions,
                                    integrationOrder=self.integrator.polynomialOrder)
        return copy

    
    def integrateStiffness( self, globalStiffnessMatrix ):
        '''
        Integrate conductivity * dNi * dNj
        '''
        basisDerivativeCache    = [ self.integrator.createCache( self.basisDerivatives[i], self.basisDerivatives.domain ) for i in range(len(self.basisDerivatives)) ]
        
        for i in range(len(self.basisDerivatives)):
            for j in range( i, len(self.basisDerivatives) ):
                value = self._invJacobian \
                        * self.integrator.integrateCached(  lambda x: self.conductivity, 
                                                            basisDerivativeCache[i]*basisDerivativeCache[j],
                                                            self.basisDerivatives.domain )
                self.updateGlobalMatrix( globalStiffnessMatrix, i, j, value )


    def integrateMass( self, globalMassMatrix ):
        '''
        Integrate capacity * Ni * Nj
        '''
        basisCache  = [ self.integrator.createCache( self.basisFunctions[i], self.basisFunctions.domain ) for i in range(len(self.basisFunctions)) ]

        for i in range(len(self.basisFunctions)):
            for j in range( i, len(self.basisFunctions) ):
                value = self._jacobian \
                        * self.integrator.integrateCached(  lambda x: self.capacity, 
                                                            basisCache[i]*basisCache[j],
                                                            self.basisFunctions.domain )
                self.updateGlobalMatrix( globalMassMatrix, i, j, value )


    def integrateLoad( self, globalLoadVector ):
        '''
        Integrate load * Ni
        '''
        cache = self.integrator.createCache( lambda x: self.load(self.toGlobalCoordinates(x)), self.basisFunctions.domain )
        for i, function in enumerate(self.basisFunctions):
            value = self._jacobian * self.integrator.integrateCached(   lambda x: function(x), 
                                                                        cache,
                                                                        self.basisFunctions.domain )
            self.updateGlobalVector( globalLoadVector, i, value )



class NonlinearHeatElement1D( Element1D ):
    '''
    Basic heat element with temperature-dependent coefficients
    '''
    def __init__( self, capacity, capacityDerivative, conductivity, conductivityDerivative, *args, **kwargs ):
        Element1D.__init__( self, *args, **kwargs )
        if not callable(capacity) or not callable(conductivity):
            raise AttributeError( "The coefficients of NonlinearHeatElement1D must be callable (functions of temperature)" )

        self.capacity               = capacity
        self.capacityDerivative     = capacityDerivative
        self.conductivity           = conductivity
        self.conductivityDerivative = conductivityDerivative


    def copy( self ):
        copy = NonlinearHeatElement1D(  self.capacity,
                                        self.capacityDerivative,
                                        self.conductivity,
                                        self.conductivityDerivative,
                                        self.domain,
                                        self.DoFs,
                                        self.load,
                                        basisFunctions=self.basisFunctions,
                                        integrationOrder=self.integrator.polynomialOrder)
        return copy


    def integrateStiffness( self, globalStiffnessMatrix, temperatureFunction, solution ):
        '''
        Integrate conductivity * dNi * dNj
        '''
        # Create caches
        if self._basisDerivativeCache is None:
            self._basisDerivativeCache    = [ self.integrator.createCache( function, self.basisDerivatives.domain ) for function in self.basisDerivatives ]
        conductivity            = lambda xi: self.conductivity( temperatureFunction(self.toGlobalCoordinates(xi)) )
        conductivityCache       = self.integrator.createCache( conductivity, self.basisFunctions.domain )
        
        # Integrate
        for i in range(len(self.basisDerivatives)):
            cache   = conductivityCache * self._basisDerivativeCache[i]
            for j in range( i, len(self.basisDerivatives) ):
                value   = self._invJacobian * self.integrator.integrateCached(  lambda x: 1.0, 
                                                                                cache * self._basisDerivativeCache[j],
                                                                                self.basisDerivatives.domain )
                self.updateGlobalMatrix( globalStiffnessMatrix, i, j, value )


    def integrateMass( self, globalMassMatrix, temperatureFunction, solution ):
        '''
        Integrate capacity * Ni * Nj
        '''
        # Create caches
        if self._basisCache is None:
            self._basisCache      = [ self.integrator.createCache( function, self.basisFunctions.domain ) for function in self.basisFunctions ]
        capacityCache   = self.integrator.createCache( lambda x: self.capacity(temperatureFunction(self.toGlobalCoordinates(x))), self.basisFunctions.domain )

        # Integrate
        for i in range(len(self.basisFunctions)):
            cache = capacityCache * self._basisCache[i]
            for j in range( i, len(self.basisFunctions) ):
                value   = self._jacobian * self.integrator.integrateCached( lambda x: 1.0,
                                                                            cache * self._basisCache[j], 
                                                                            self.basisFunctions.domain )
                self.updateGlobalMatrix( globalMassMatrix, i, j, value )


    def integrateLoad( self, globalLoadVector, temperatureFunction, solution ):
        '''
        Integrate load * Ni
        '''
        # Create caches
        if self._basisCache is None:
            self._basisCache              = [ self.integrator.createCache( function, self.basisFunctions.domain ) for function in self.basisFunctions ]
        loadCache = self.integrator.createCache( lambda x: self.load(self.toGlobalCoordinates(x)), self.basisFunctions.domain )
        
        # Integrate
        for i in range(len(self.basisFunctions)):
            value   = self._jacobian * self.integrator.integrateCached( lambda x: 1.0, 
                                                                        loadCache * self._basisCache[i],
                                                                        self.basisFunctions.domain )
            self.updateGlobalVector( globalLoadVector, i, value )


    def integrateGeometricStiffness( self, matrix, temperatureFunction, solution ):
        '''
        Integrate conductivityDerivative * dNi*Nj * sum(dNk*uk)
        '''
        # Create caches
        if self._basisDerivativeCache is None:
            self._basisDerivativeCache    = [ self.integrator.createCache( function, self.basisDerivatives.domain ) for function in self.basisDerivatives ]
        if self._basisCache is None:
            self._basisCache              = [ self.integrator.createCache( function, self.basisFunctions.domain ) for function in self.basisFunctions ]
        dConductivity           = lambda xi: self.conductivityDerivative( temperatureFunction(self.toGlobalCoordinates(xi)) )
        dConductivityCache      = self.integrator.createCache( dConductivity, self.basisDerivatives.domain )
        
        # Integrate
        for i in range(len( self.basisDerivatives )):
            cache0 = self._basisDerivativeCache[i] * dConductivityCache
            for j in range(len( self.basisFunctions )):
                cache1  = cache0 * self._basisCache[j]
                value   = 0.0
                for k in range(len( self.basisDerivatives )):
                    incr    =  self._invJacobian * self.integrator.integrateCached(    lambda x: 1.0,
                                                                                        cache1 * self._basisDerivativeCache[k],
                                                                                        self.basisFunctions.domain )
                    value   += solution[self.DoFs[k]] * incr

                matrix[ self.DoFs[i], self.DoFs[j] ] += value