# --- Python Imports ---
import numpy as np
import scipy.sparse as sparse
import functools
from copy import copy
from enum import Enum

# --- Internal Imports ---
from cie.fem.utilities import isNumpyArray, gridPairs
from .element import Element, Element1D

# ---------------------------------------------------------
# Decorators
def requiresInitialized( function ):
    @functools.wraps( function )
    def decoratedFunction( instance, *args, **kwargs ):
        for key, value in instance.__dict__.items():
            if value is None:
                raise AttributeError( "Uninitialized member variable when calling a function that requires an initialized object!:\n" + key )
        else:
            return function( instance, *args, **kwargs )
    return decoratedFunction

# ---------------------------------------------------------
class FEModel:
    '''
    Aggregate base class to hold elements, store global matrices, and apply constraints to
    '''
    def __init__( self, size ):
        self.size   = size
        self.shape  = np.asarray( (size, size), dtype=np.uint64 )

        # Initialize elements
        self.elements       = None

        # Initialize matrices
        self.stiffness      = None
        self.load           = None

        self.boundaries     = []


    def resetMatrices( self, stiffness=True, load=True ):
        # Check arguments: stiffness and load must be both True or both False
        #if stiffness is not load:
        #    raise RuntimeError( "The stiffness matrix and load vector must be reset together, \
        #                        otherwise boundary conditions won't be enforced correctly" )

        # Reset structural matrices
        if stiffness:
            self.stiffness.data = np.zeros( self.stiffness.data.shape, dtype=self.stiffness.data.dtype )

        if load:
            self.load           = np.zeros( self.load.shape, dtype=self.load.dtype )

        # Reset boundaries
        if stiffness or load:
            for boundary in self.boundaries:
                boundary.applied = False


    def getDoFs( self ):
        '''
        Get a unique set of DoF pairs that indicate nonzero entries in matrices
        '''
        if self.elements is None:
            raise AttributeError( "Attempt to get DoFs of an uninitialized model" )

        DoFs = None
        for element in self.elements:
            if DoFs is None:
                DoFs = gridPairs( element.DoFs )
            else:
                DoFs = np.concatenate( (DoFs, gridPairs( element.DoFs )), axis=0 )

        return np.unique( DoFs, axis=0 )


    def allocateZeros( self ):
        # Get nonzero indices
        DoFs = np.transpose( self.getDoFs(), (1,0) )

        # Check if already initialized
        if self.stiffness is not None:
            raise RuntimeError( "Attempt to initialize existing matrices" )

        # Initialize matrices
        self.stiffness  = sparse.csr_matrix(    ( np.zeros( DoFs.shape[1] ), DoFs ),
                                                shape=self.shape )
        self.load       = np.zeros( self.size )

        return DoFs


    @requiresInitialized
    def integrate( self, stiffness=True, load=True ):
        '''
        Loop through elements and integrate their matrices
        '''
        for element in self.elements:
            if stiffness:
                element.integrateStiffness( self.stiffness )
            if load:
                element.integrateLoad( self.load )


    def addBoundaryCondition( self, boundaryCondition, *args, **kwargs ):
        self.boundaries.append( boundaryCondition )
        self.applyBoundaryCondition( boundaryCondition, *args, **kwargs )
        return len(self.boundaries)-1


    @requiresInitialized
    def applyBoundaryCondition( self, boundaryCondition, *args, **kwargs ):
        # Check if boundary condition is initialized
        boundaryCondition.initialize(self)

        # Defer enforcing to dedicated functions
        if boundaryCondition.BCType == "dirichlet":
            self.applyDirichletBoundary( boundaryCondition, *args, **kwargs )
        elif boundaryCondition.BCType == "neumann":
            self.applyNeumannBoundary( boundaryCondition, *args, **kwargs )
        else:
            raise NotImplementedError( "Unknown boundary condition type" )


    @requiresInitialized
    def applyDirichletBoundary( self, boundaryCondition ):
        if boundaryCondition.applied is not True:
            # Reset the affected DoF
            for index in range(self.size):
                if np.abs( self.stiffness[boundaryCondition.DoF, index] )>1e-15:
                    self.stiffness[boundaryCondition.DoF, index] = 0.0
                if np.abs( self.stiffness[ index, boundaryCondition.DoF] )>1e-15:
                    self.load[index] -= self.stiffness[index, boundaryCondition.DoF] * boundaryCondition.value
                    self.stiffness[ index, boundaryCondition.DoF ] = 0.0

            # Enforce boundary
            self.stiffness[boundaryCondition.DoF,boundaryCondition.DoF] = 1.0
            self.load[boundaryCondition.DoF] = boundaryCondition.value
            #self.stiffness[boundaryCondition.DoF,boundaryCondition.DoF] += boundaryCondition.penaltyValue
            #self.load[boundaryCondition.DoF]    += boundaryCondition.penaltyValue * boundaryCondition.value

            boundaryCondition.applied = True


    @requiresInitialized
    def applyNeumannBoundary( self, boundaryCondition ):
        self.load[boundaryCondition.DoF]    = boundaryCondition.value
        boundaryCondition.applied           = True


    @requiresInitialized
    def sample( self, solution, samples ):
        values = np.zeros( samples.shape )
        for element in self.elements:
            values += element( solution[element.DoFs], samples )

        return values


    @requiresInitialized
    def sampleDerivative( self, solution, samples ):
        values = np.zeros( samples.shape )
        for element in self.elements:
            values += element.derivative( solution[element.DoFs], samples )

        return values


    @requiresInitialized
    def updateLoad( self, load, *args, **kwargs ):
        '''
        Set load on all elements and reintegrate
        '''
        # Clear
        self.load = np.zeros( self.load.shape )

        # Set load on all elements and integrate
        for element in self.elements:
            element.load = load
            element.integrateLoad( self.load, *args, **kwargs )




class TransientFEModel( FEModel ):
    '''
    FEModel extended with mass matrix and time-and-space-dependent load function
    '''
    def __init__( self, *args, loadFunction=None, **kwargs ):
        FEModel.__init__( self, *args, **kwargs )

        # Initialize
        self.mass           = None
        self._loadFunction  = loadFunction
        self.time           = 0.0


    def resetMatrices( self, mass=True, **kwargs ):
        FEModel.resetMatrices(self, **kwargs)
        if mass:
            self.mass.data      = np.zeros( self.mass.data.shape, dtype=self.mass.data.dtype )


    @property
    def loadFunction( self ):
        if self._loadFunction is None:
            raise AttributeError( "Unset load function!" )

        return self._loadFunction


    def allocateZeros( self ):
        DoFs = FEModel.allocateZeros( self )

        if self.mass is not None:
            raise RuntimeError( "Attempt to initialize existing matrices" )
        self.mass   = sparse.csc_matrix((   np.zeros( DoFs.shape[1] ), DoFs ),
                                            shape=self.shape )
        return DoFs


    @requiresInitialized
    def integrate( self, *args, stiffness=True, mass=True, load=True, **kwargs ):
        FEModel.integrate( self, *args, stiffness=stiffness, load=load, **kwargs )
        if mass:
            for element in self.elements:
                    element.integrateMass( self.mass, *args, **kwargs )


    @requiresInitialized
    def updateTime( self, time ):
        self.time = time
        self.updateLoad( lambda x: self.loadFunction( time, x ) )
        for boundaryCondition in self.boundaries:
            self.applyBoundaryCondition( boundaryCondition )


    @requiresInitialized
    def addBoundaryCondition( self, boundaryCondition ):
        if callable(boundaryCondition.value) is False:
            raise AttributeError( "Boundary condition value for a TransientFEModel must be a function of time!" )
        return FEModel.addBoundaryCondition( self, boundaryCondition )


    @requiresInitialized
    def applyBoundaryCondition( self, boundaryCondition ):
        BC                          = copy(boundaryCondition)
        BC.value                    = BC.value( self.time )
        boundaryCondition.applied   = False
        return FEModel.applyBoundaryCondition( self, BC )


    @requiresInitialized
    def removeBoundaryCondition( self, boundaryCondition ):
        BC                          = copy(boundaryCondition)
        BC.value                    = BC.value( self.time )
        boundaryCondition.applied   = False
        return FEModel.removeBoundaryCondition( self, BC )




class NonlinearFEModel( FEModel ):
    '''
    FEModel with material nonlinearities
    '''
    def __init__( self, *args, **kwargs ):
        FEModel.__init__( self, *args, **kwargs )
        self.geometricStiffness = None


    def resetMatrices( self, geometricStiffness=True, **kwargs ):
        FEModel.resetMatrices(self, **kwargs)
        if geometricStiffness:
            self.geometricStiffness.data = np.zeros(    self.geometricStiffness.data.shape,
                                                        self.geometricStiffness.data.dtype )


    def allocateZeros( self ):
        DoFs                    = FEModel.allocateZeros( self )
        self.geometricStiffness = self.stiffness.copy()
        return DoFs


    @requiresInitialized
    def integrate( self, *args, stiffness=True, load=True, geometricStiffness=True, **kwargs ):
        for element in self.elements:
            if stiffness:
                element.integrateStiffness( self.stiffness, *args, **kwargs )
            if load:
                element.integrateLoad( self.load, *args, **kwargs )
            if geometricStiffness:
                element.integrateGeometricStiffness( self.geometricStiffness, *args, **kwargs )


    @requiresInitialized
    def applyDirichletBoundary(self, boundaryCondition):
        if not boundaryCondition.applied:
            FEModel.applyDirichletBoundary(self, boundaryCondition)
            for i in range(self.size):
                if np.abs(self.geometricStiffness[boundaryCondition.DoF, i]) > 1e-15:
                    self.geometricStiffness[boundaryCondition.DoF, i] = 0.0



    #@requiresInitialized
    #def applyNeumannBoundary( self, boundaryCondition, fieldValue ):
    #    # TODO
    #    FEModel.applyNeumannBoundary( self, boundaryCondition )



class NonlinearTransientFEModel( TransientFEModel ):
    def __init__( self, *args, **kwargs ):
        TransientFEModel.__init__( self, *args, **kwargs )
        self.geometricStiffness = None


    def resetMatrices( self, geometricStiffness=True, **kwargs ):
        TransientFEModel.resetMatrices(self, **kwargs)
        if geometricStiffness:
            self.geometricStiffness.data = np.zeros(    self.geometricStiffness.data.shape,
                                                        self.geometricStiffness.data.dtype )


    def allocateZeros( self ):
        DoFs                    = TransientFEModel.allocateZeros( self )
        self.geometricStiffness = self.stiffness.copy()
        return DoFs


    @requiresInitialized
    def integrate( self, *args, stiffness=True, mass=True, load=True, geometricStiffness=True, **kwargs ):
        for element in self.elements:
            if stiffness:
                element.integrateStiffness( self.stiffness, *args, **kwargs )
            if load:
                element.integrateLoad( self.load, *args, **kwargs )
            if mass:
                element.integrateMass( self.mass, *args, **kwargs )
            if geometricStiffness:
                element.integrateGeometricStiffness( self.geometricStiffness, *args, **kwargs )


    @requiresInitialized
    def updateTime( self, time, solution, stiffness=True, mass=True, geometricStiffness=True ):
        '''
        Update model time
        Triggers the following:
            - the load function gets updated to the current time
            - boundary conditions get updated to the current time
            - reintegrates the load function
            - upon request, reintegrates the stiffness, mass, and geometric stiffness matrices
        '''
        # Update time
        self.time = time

        # Update load function and reintegrate load vector
        self.updateLoad(    lambda x: self.loadFunction(time, x),
                            lambda x: self.sample( solution, x ),
                            solution )

        # Reintegrate structural matrices if requested
        self.resetMatrices( stiffness=stiffness,
                            mass=mass,
                            geometricStiffness=geometricStiffness,
                            load=False )
        self.integrate( lambda x: self.sample( solution, x ),
                        solution,
                        stiffness=stiffness,
                        mass=mass,
                        geometricStiffness=geometricStiffness,
                        load=False )

        # Apply boundary conditions
        for boundaryCondition in self.boundaries:
            self.applyBoundaryCondition( boundaryCondition )


    @requiresInitialized
    def applyDirichletBoundary(self, boundaryCondition):
        if not boundaryCondition.applied:
            FEModel.applyDirichletBoundary(self, boundaryCondition)
            for i in range(self.size):
                if np.abs(self.geometricStiffness[boundaryCondition.DoF, i]) > 1e-15:
                    self.geometricStiffness[boundaryCondition.DoF, i] = 0.0