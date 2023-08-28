
# ---------------------------------------------------------
class BoundaryCondition:
    '''
    Boundary condition base class
    '''
    def __init__( self, DoF, position, value ):
        self.DoF            = DoF
        self.position       = position
        self.value          = value
        
        self.elementIDs     = None          #<-- Indices of elements in the FEModel that affect the DoF of the boundary condition
        self.BCType         = None
        self.applied        = False


    def initialize(self, model):
        '''
        Fill necessary members
        '''
        if self.elementIDs is None:
            self.elementIDs = []
            for elementID, element in enumerate(model.elements):
                if self.DoF in element.DoFs:
                    self.elementIDs.append(elementID)


class PenaltyBoundaryCondition( BoundaryCondition ):
    def __init__( self, *args, penaltyValue=1e5, **kwargs ):
        BoundaryCondition.__init__( self, *args, **kwargs )
        self.penaltyValue   = penaltyValue


class DirichletBoundary( PenaltyBoundaryCondition ):
    def __init__( self, *args, **kwargs ):
        PenaltyBoundaryCondition.__init__( self, *args, **kwargs )
        self.BCType = "dirichlet"


class NeumannBoundary( BoundaryCondition ):
    def __init__( self, *args, **kwargs ):
        BoundaryCondition.__init__( self, *args, **kwargs )
        self.BCType = "neumann"