# --- Python Imports ---
import numpy as np
import scipy.sparse.linalg as linalg
from matplotlib import pyplot as plt

# --- Internal Imports ---
from cie.fem.discretization import IntegratedHierarchicBasisFunctions
from cie.fem.discretization import LinearHeatElement1D
from cie.fem.discretization import FEModel
from cie.fem.discretization import DirichletBoundary, NeumannBoundary
from cie.fem.numeric import solveLinearSystem

# ---------------------------------------------------------
# Geometry and material
length              = 1.0
capacity            = 1.0
conductivity        = 1.0

# Load
#load                = lambda x: np.sin(x*np.pi)
load                = lambda x: 0.0

# Discretization
nElements           = 5
polynomialOrder     = 1

# Integration
integrationOrder    = 2*polynomialOrder + 1

# ---------------------------------------------------------
# Initialize FE model
model               = FEModel( nElements*polynomialOrder + 1 )

# Create elements
basisFunctions      = IntegratedHierarchicBasisFunctions( polynomialOrder=polynomialOrder )
model.elements      = [ LinearHeatElement1D(    capacity,
                                                conductivity,
                                                (i*length/nElements, (i+1)*length/nElements),
                                                np.asarray( range(i*polynomialOrder, (i+1)*polynomialOrder+1) ),
                                                load,
                                                basisFunctions=basisFunctions,
                                                integrationOrder=integrationOrder   ) 
                        for i in range(nElements) ]

# Integrate
model.allocateZeros( )
model.integrate( )

# Boundary conditions
leftBCID    = model.addBoundaryCondition(   DirichletBoundary(  0, 
                                                                0.0,
                                                                1.0 ))

rightBCID   = model.addBoundaryCondition(   NeumannBoundary(  nElements*polynomialOrder,
                                                                length,
                                                                1.0) )

# Solve
u       = solveLinearSystem( model.stiffness, model.load )

# Output
np.set_printoptions( precision=1, suppress=True )
print( model.stiffness.todense() )
print( "Solution: " + str(u) )

samples = np.linspace( 0, length, num=100 )
values  = model.sample( u, samples )
plt.plot( samples, values, '.-' )
plt.show()