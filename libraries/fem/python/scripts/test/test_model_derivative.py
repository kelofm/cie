# --- Python Imports ---
import numpy as np
import scipy.sparse.linalg as linalg
import scipy
from matplotlib import pyplot as plt

# --- Internal Imports ---
from cie.fem.discretization import IntegratedHierarchicBasisFunctions
from cie.fem.discretization import NonlinearHeatElement1D
from cie.fem.discretization import NonlinearFEModel
from cie.fem.discretization import DirichletBoundary, NeumannBoundary
from cie.fem.numeric import stationaryLoadControl as nonlinearSolver
from cie.fem.postprocessing import ConvergencePlot

# ---------------------------------------------------------
# Geometry and material
length              = 1.0
capacity            = lambda u: 1.0
dCapacity           = lambda u: 0.0
#conductivity        = lambda u: 1.0
#dConductivity       = lambda u: 0.0
conductivity        = lambda u: 1.0 + 9.0 * np.exp( -(u-0.5)**2 / 0.005 )
dConductivity       = lambda u: 9.0 * np.exp( -(u-0.5)**2 / 0.005 ) * (2.0/0.005)*(0.5-u)

# Load
load                = lambda x: 0.0
boundaryFlux        = 2.0

# Discretization
nElements           = 10
polynomialOrder     = 3

# Integration
integrationOrder    = 3 * (2*polynomialOrder + 1)

# Iteration
baseIncrement       = 0.2
minIncrement        = 0.01
maxIncrement        = 0.4
maxIncrements       = 15
maxCorrections      = 8
tolerance           = 1e-5

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

# Boundary conditions (right BC is the load)
leftBCID    = model.addBoundaryCondition(   DirichletBoundary(  0, 
                                                                0.0,
                                                                0.0 ))

rightBCID   = model.addBoundaryCondition(   NeumannBoundary(    nElements*polynomialOrder,
                                                                length,
                                                                boundaryFlux) )

# Solve
u = nonlinearSolver(    model,
                        np.zeros(model.size),
                        baseIncrement=baseIncrement,
                        minIncrement=minIncrement,
                        maxIncrement=maxIncrement,
                        maxCorrections=maxCorrections,
                        tolerance=tolerance,
                        verbose=False   )

# ---------------------------------------------------------
# Approximate derivative
dx                      = 1e-12
approximateDerivative   = lambda x: ( model.sample(u,x+dx) - model.sample(u,x) ) / dx
derivative              = lambda x: model.sampleDerivative(u,x)
squaredDifference       = lambda x: (approximateDerivative(x)-derivative(x))**2

samples                 = np.linspace(0.0,length-dx,num=100)
err                     = np.trapz( [squaredDifference(x) for x in samples],
                                    x=samples )

print( "Error: %.5f" % err )

#if err > 1e-3:
#    plt.plot( samples, model.sample(u,samples), "k.-" )
#    plt.plot( samples, [approximateDerivative(x) for x in samples], "b" )
#    plt.plot( samples, [derivative(x) for x in samples], "r" )
#    plt.show()

assert( err < 1e-3 )