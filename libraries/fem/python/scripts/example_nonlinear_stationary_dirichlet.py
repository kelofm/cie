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
#conductivity        = lambda u: 1.0
#dConductivity       = lambda u: 0.0
conductivity        = lambda u: 1.0 + 9.0 * np.exp( (-(u-0.5)**2) / 0.005 )
dConductivity       = lambda u: -9.0 * np.exp( (-(u-0.5)**2) / 0.005 ) * (2.0/0.005)*(u-0.5)

# Load
load                = lambda x: 0.0
boundaryTemperature = 1.0

# Discretization
nElements           = 10
polynomialOrder     = 3

# Integration
integrationOrder    = 2 * (2*polynomialOrder + 1)

# Iteration
baseIncrement       = 0.2
minIncrement        = 0.01
maxIncrement        = 0.4
maxIncrements       = 15
maxCorrections      = 8
tolerance           = 1e-5

# ---------------------------------------------------------
# General initialization
samples         = np.linspace( 0, length, num=100 )
fig             = plt.figure( )
axes            = ( fig.add_subplot( 2,1,1 ),
                    fig.add_subplot( 2,1,2 ))

convergencePlot = ConvergencePlot()

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

rightBCID   = model.addBoundaryCondition(   DirichletBoundary(  nElements*polynomialOrder,
                                                                length,
                                                                boundaryTemperature) )

# ---------------------------------------------------------
# Solve
u = nonlinearSolver(    model,
                        np.zeros(model.size),
                        baseIncrement=baseIncrement,
                        minIncrement=minIncrement,
                        maxIncrement=maxIncrement,
                        maxIncrements=maxIncrements,
                        maxCorrections=maxCorrections,
                        tolerance=tolerance,
                        verbose=True,
                        axes=axes[1],
                        convergencePlot=convergencePlot   )

# ---------------------------------------------------------
# Plot conductivity
samples = np.linspace( 0.0, np.max(model.sample( u, samples )), num=len(samples) )
axes[0].plot( samples, [conductivity(temp) for temp in samples] )
axes[0].set_xlabel( "T [C]" )
axes[0].set_ylabel( r'$\kappa$' + " [W/K]" )
axes[0].set_title( "Conduction(temperature)" )

plt.tight_layout()
plt.show()