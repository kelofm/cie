# --- Python Imports ---
import numpy as np
import scipy.sparse.linalg as linalg
from matplotlib import pyplot as plt

# --- Internal Imports ---
from cie.fem.discretization import IntegratedHierarchicBasisFunctions
from cie.fem.discretization import NonlinearHeatElement1D
from cie.fem.discretization import NonlinearFEModel
from cie.fem.discretization import DirichletBoundary
from cie.fem.numeric import stationaryLoadControl as nonlinearSolver
from cie.fem.postprocessing import ConvergencePlot

# ---------------------------------------------------------
# Geometry and material
length              = 1.0
capacity            = lambda u: 1.0
dCapacity           = lambda u: 0.0
#conductivity        = lambda u: 1.0 + u*u
#dConductivity       = lambda u: 2.0*u
conductivity        = lambda u: 1.0 + 9.0 * np.exp( -(u-0.5)**2 / 0.005 )
dConductivity       = lambda u: 9.0 * np.exp( -(u-0.5)**2 / 0.005 ) * (2.0/0.005)*(0.5-u)


# Load
load                = lambda x: 0.0
boundaryTemperature = 1.0

# Discretization
nElements           = 5
polynomialOrder     = 1

# Integration
integrationOrder    = 2 * (2*polynomialOrder + 1)

# ---------------------------------------------------------
# General initialization
samples         = np.linspace( 0, length, num=100 )

# ---------------------------------------------------------
# Initialize FE model
model               = NonlinearFEModel( nElements*polynomialOrder + 1 )

# Create elements
basisFunctions      = IntegratedHierarchicBasisFunctions( polynomialOrder=polynomialOrder )
model.elements      = [ NonlinearHeatElement1D( capacity,
                                                dCapacity,
                                                conductivity,
                                                dConductivity,
                                                (i*length/nElements, (i+1)*length/nElements),
                                                np.asarray( range(i*polynomialOrder, (i+1)*polynomialOrder+1) ),
                                                load,
                                                basisFunctions=basisFunctions,
                                                integrationOrder=integrationOrder   ) 
                        for i in range(nElements) ]

model.allocateZeros( )
derivative  = model.stiffness.copy()

# Boundary conditions (right BC is the load)
leftBCID    = model.addBoundaryCondition(   DirichletBoundary(  0, 
                                                                0.0,
                                                                0.0 ))

rightBCID   = model.addBoundaryCondition(   DirichletBoundary(  nElements*polynomialOrder,
                                                                length,
                                                                boundaryTemperature) )

# ---------------------------------------------------------
u           = np.random.rand(model.size)
u[0]        = 0.0
u[-1]       = boundaryTemperature
du          = 1e-6
dr          = []

# Compute analytical derivative
model.resetMatrices()
model.integrate(    lambda x: model.sample(u, x),
                    u )
for boundary in model.boundaries:
    model.applyBoundaryCondition(boundary)

tangentStiffness    = model.stiffness + model.geometricStiffness
r                   = model.stiffness.dot(u) - model.load

# Compute finite difference derivative
for i in range(model.size):
    u[i]    += du

    # Recompute structural matrices
    model.resetMatrices()
    model.integrate(    lambda x: model.sample(u, x),
                        u )
    for boundary in model.boundaries:
        model.applyBoundaryCondition(boundary)

    dr.append( (model.stiffness.dot(u) - model.load - r)/du )
    
    u[i]    -= du

# ---------------------------------------------------------
dr  = np.transpose( np.asarray(dr), (1,0) )
np.set_printoptions( precision=2, suppress=True )

print( dr )
print( tangentStiffness.todense() )
print( tangentStiffness.todense() - dr )
