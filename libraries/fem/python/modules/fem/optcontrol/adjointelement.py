# --- Discretization Imports ---
from cie.fem.discretization.element import Element1D

# --- Python Imports ---
import numpy as np

# ---------------------------------------------------------

class AdjointHeatElement1D( Element1D ):
    '''
    Element type for AdjointModel
    The coefficients (capacity, conductivity) don't depend on the solution directly,
    but on the position, so integration must be carried out accordingly.

    The constructor expects a NonlinearHeatElement1D and a forward solution, from which
    the coefficients are constructed. An extra integration target (nonSymmetricStffness)
    is added.
    '''
    def __init__(   self, 
                    nonlinearHeatElement1D, 
                    forwardModel, 
                    forwardSolution, 
                    referenceSolution,
                    ddConductivity      # Second derivative of the original element's conductivity
                                        # TODO: create a wrapper function for these coefficients that
                                        #       encapsulate their derivatives to avoid this clutter
                     ):
        '''
        nonlinearHeatElement1D  : element from the forward model
        model                   : forward model of type NonlinearTransientFEModel
        forwardSolution         : list of coefficients representing the temperature field u(t,x)
                                    for the specified time points
        '''
        super( AdjointHeatElement1D, self ).__init__(   nonlinearHeatElement1D.domain,
                                                        nonlinearHeatElement1D.DoFs,
                                                        nonlinearHeatElement1D.load,
                                                        basisFunctions=nonlinearHeatElement1D.basisFunctions,
                                                        integrationOrder=nonlinearHeatElement1D.integrator.polynomialOrder )
        self.element            = nonlinearHeatElement1D
        self.solutionField      = forwardModel.sample
        self.derivativeField    = forwardModel.sampleDerivative
        self.forwardSolution    = forwardSolution

        self.timeIndex          = 0

        # Second derivative of the original element's conductivity
        # TODO: create a wrapper class for these coefficients that
        #       encapsulates their derivatives
        self._ddConductivity    = ddConductivity

        # TODO: Caching


    def capacity( self, x ):
        u = self.solutionField( self.forwardSolution[self.timeIndex], x )
        return self.element.capacity(u) + self.element.capacityDerivative(u)*u


    def conductivity( self, x ):
        u = self.solutionField( self.forwardSolution[self.timeIndex], x )
        return self.element.conductivity(u) + self.element.conductivityDerivative(u)*u


    def conductivityDerivative( self, x ):
        u   = self.solutionField( self.forwardSolution[self.timeIndex], x )
        du  = self.derivativeField( self.forwardSolution[self.timeIndex], x )
        return 2*self.element.conductivityDerivative(u)*du + self._ddConductivity(u)*u


    def integrateNonsymmetricStiffness( self, globalMatrix ):
        '''
        Integrate conductivityDerivative * Ni*dNj
        '''
        basisCache                  = [ self.integrator.createCache( self.basisFunctions[i], self.basisFunctions.domain ) for i in range(len(self.basisFunctions)) ]
        basisDerivativeCache        = [ self.integrator.createCache( self.basisDerivatives[i], self.basisDerivatives.domain ) for i in range(len(self.basisDerivatives)) ]
        conductivityDerivativeCache = self.integrator.createCache( lambda xi: self.conductivityDerivative(self.toGlobalCoordinates(xi)), self.basisFunctions.domain )

        for i in range(len(self.basisFunctions)):
            cache = conductivityDerivativeCache * basisCache[i]
            for j in range(len(self.basisDerivatives)):
                value = self.integrator.integrateCached(    lambda x: 1.0,
                                                            cache * basisDerivativeCache[j],
                                                            self.basisFunctions.domain )
                globalMatrix[self.DoFs[i],self.DoFs[j]] += value


    def integrateMass( self, globalMatrix ):
        '''
        Integrate capacity * Ni * Nj
        '''
        basisCache  = [ self.integrator.createCache( self.basisFunctions[i], self.basisFunctions.domain ) for i in range(len(self.basisFunctions)) ]

        for i in range(len(self.basisFunctions)):
            for j in range( i, len(self.basisFunctions) ):
                value = self._jacobian \
                        * self.integrator.integrateCached(  lambda xi: self.capacity( self.toGlobalCoordinates(xi) ), 
                                                            basisCache[i]*basisCache[j],
                                                            self.basisFunctions.domain )
                self.updateGlobalMatrix( globalMatrix, i, j, value )


    def integrateStiffness( self, globalMatrix ):
        '''
        Integrate conductivity * dNi * dNj
        '''
        basisDerivativeCache    = [ self.integrator.createCache( self.basisDerivatives[i], self.basisDerivatives.domain ) for i in range(len(self.basisDerivatives)) ]
        
        for i in range(len(self.basisDerivatives)):
            for j in range( i, len(self.basisDerivatives) ):
                value = self._invJacobian \
                        * self.integrator.integrateCached(  lambda xi: self.conductivity( self.toGlobalCoordinates(xi) ), 
                                                            basisDerivativeCache[i]*basisDerivativeCache[j],
                                                            self.basisDerivatives.domain )
                self.updateGlobalMatrix( globalMatrix, i, j, value )


    def integrateLoad( self, globalVector ):
        '''
        Not element-dependent
        '''
        raise RuntimeError( "The adjoint load vector is not computed element-wise" )