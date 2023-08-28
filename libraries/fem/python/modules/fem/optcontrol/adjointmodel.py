# --- Discretization Imports ---
from cie.fem.discretization.element import Element1D
from cie.fem.discretization import TransientFEModel, NonlinearTransientFEModel

# --- Internal Imports ---
from .adjointelement import AdjointHeatElement1D

# --- Python Imports ---
import numpy as np
import scipy.sparse as sparse

# ---------------------------------------------------------



class AdjointModel(TransientFEModel):
    '''
    Nonlinear adjoint model attached to a NonlinearTransientFEModel
    '''
    def __init__( self, model, time, referenceTimeSeries, forwardTimeSeries, ddConductivity, *args, **kwargs ):
        TransientFEModel.__init__( self, model.size )

        # Check model type
        # The adjoint system is linear with respect to its state, the system
        # it's derived from must be nonlinear. A simpler structure is used
        # for solving linear systems.
        assert type(model) is NonlinearTransientFEModel, "Invalid model type!"

        self.femodel                = model
        self.time                   = time
        self.referenceSolution      = referenceTimeSeries
        self.forwardSolution        = forwardTimeSeries

        self.timeIndex              = 0
        self._load                  = None

        self.nonsymmetricStiffness  = None

        # Transform model members 
        # TODO: ddConductivity
        self.elements               = []
        for element in model.elements:
            self.elements.append(
                AdjointHeatElement1D(   element,
                                        self.femodel,
                                        self.forwardSolution,
                                        self.referenceSolution,
                                        ddConductivity ) )

        self.computeLoad()


    def allocateZeros( self ):
        DoFs = TransientFEModel.allocateZeros(self)
        self.nonsymmetricStiffness  = sparse.csr_matrix(    ( np.zeros( DoFs.shape[1] ), DoFs ), 
                                                            shape=self.shape )


    def resetMatrices( self, stiffness=True, mass=True, load=False ):
        TransientFEModel.resetMatrices( self, stiffness=stiffness, mass=mass, load=load )
        if stiffness:
            self.nonsymmetricStiffness.data = np.zeros( self.nonsymmetricStiffness.data.shape, dtype=self.nonsymmetricStiffness.data.dtype )


    def computeLoad( self ):
        '''
        Integrate Ni*Nj*(uref - uforward)
        '''
        DoFs = np.transpose( self.getDoFs(), (1,0) )
        bareMass    = sparse.csr_matrix(    ( np.zeros( DoFs.shape[1] ), DoFs ), 
                                            shape=self.shape )
        for element in self.elements:
            basisCache = [ element.integrator.createCache( element.basisFunctions[i], element.basisFunctions.domain ) for i in range(len(element.basisFunctions)) ]
            for i in range(len(element.basisFunctions)):
                for j in range( i, len(element.basisFunctions) ):
                    value = element._jacobian \
                            * element.integrator.integrateCached(   lambda x: 1.0, 
                                                                    basisCache[i]*basisCache[j],
                                                                    element.basisFunctions.domain )
                    element.updateGlobalMatrix( bareMass, i, j, value )

        self._load = [ bareMass.dot( reference - forward ) for reference,forward in zip(self.referenceSolution,self.forwardSolution) ]
        


    @property
    def load( self ):
        '''
        Adjoint load vector: depends only on the reference and forward solutions,
        and doesn't have to be integrated
        '''
        return self._load[ self.timeIndex ]

    
    @load.setter
    def load( self, value ):
        pass


    @property
    def _loadFunction( self ):
        raise RuntimeError( "_loadFunction and loadFunction should not be used for AdjointModel" )


    @_loadFunction.setter
    def _loadFunction( self, value ):
        pass


    def updateTime( self, timeIndex ):
        self.timeIndex = timeIndex
        for element in self.elements:
            element.timeIndex = timeIndex
        self.integrate()
        for boundaryCondition in self.boundaries:
            self.applyBoundaryCondition(boundaryCondition)


    def integrate( self, *args, stiffness=True, mass=True, load=False, **kwargs ):
        self.resetMatrices( stiffness=stiffness, mass=mass, load=load )
        TransientFEModel.integrate(self,*args,stiffness=stiffness,mass=mass,load=load,**kwargs)
        if stiffness:
            for element in self.elements:
                element.integrateNonsymmetricStiffness(self.nonsymmetricStiffness)


    