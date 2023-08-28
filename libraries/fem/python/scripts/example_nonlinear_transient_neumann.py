# --- Python Imports ---
import numpy as np
import scipy.sparse.linalg as linalg
from matplotlib import pyplot as plt

# --- Internal Imports ---
from cie.fem.discretization import IntegratedHierarchicBasisFunctions
from cie.fem.discretization import NonlinearHeatElement1D
from cie.fem.discretization import NonlinearTransientFEModel as Model
from cie.fem.discretization import DirichletBoundary, NeumannBoundary
from cie.fem.numeric import solveNonlinearHeat1D
from cie.fem.numeric import transientLoadControl as nonlinearSolver
from cie.fem.postprocessing import ConvergencePlot
from cie.fem.postprocessing.graphics import animateTimeSeries

# ---------------------------------------------------------
# Geometry and material
length                      = 1.0
#capacity                    = lambda u: 1.0
#dCapacity                   = lambda u: 0.0
capacity                    = lambda u: 1.0 + 9.0 * np.exp( -(u-0.5)**2 / 0.005 )
dCapacity                   = lambda u: 9.0 * np.exp( -(u-0.5)**2 / 0.005 ) * (2.0/0.005)*(0.5-u)
#conductivity                = lambda u: 1.0
#dConductivity               = lambda u: 0.0
conductivity                = lambda u: 1.0 + 9.0 * np.exp( -(u-0.5)**2 / 0.005 )
dConductivity               = lambda u: 9.0 * np.exp( -(u-0.5)**2 / 0.005 ) * (2.0/0.005)*(0.5-u)

# Load
load                        = lambda t, x: 0.0
boundaryFlux                = lambda t: 1.0 - np.exp(-100.0*t)
penaltyValue                = 1e10

# Discretization
time                        = np.linspace(0.0, 1.0, num=30)
nElements                   = 10
polynomialOrder             = 3

# Integration
integrationOrder            = 1 * (2*polynomialOrder + 1)
finiteDifferenceImplicity   = 0.5

# Iteration
baseIncrement       = 1.0
minIncrement        = 0.01
maxIncrement        = 1.0
maxIncrements       = 100
maxCorrections      = 15
tolerance           = 1e-10

# ---------------------------------------------------------
# General initialization
samples         = np.linspace( 0, length, num=100 )
convergencePlot = ConvergencePlot()

# ---------------------------------------------------------
# Initialize FE model
model               = Model( nElements*polynomialOrder + 1, loadFunction=load )

# Create elements
basisFunctions      = IntegratedHierarchicBasisFunctions( polynomialOrder=polynomialOrder )
initialLoad         = lambda x: load( time[0], x )
model.elements      = [ NonlinearHeatElement1D( capacity, dCapacity,
                                                conductivity, dConductivity,
                                                (i*length/nElements, (i+1)*length/nElements),
                                                np.asarray( range(i*polynomialOrder, (i+1)*polynomialOrder+1) ),
                                                initialLoad,
                                                basisFunctions=basisFunctions,
                                                integrationOrder=integrationOrder   ) 
                        for i in range(nElements) ]

# Integrate
model.allocateZeros( )

# Boundary conditions
model.addBoundaryCondition( DirichletBoundary(  0, 
                                                0.0,
                                                lambda t: 0.0   ) )

model.addBoundaryCondition( NeumannBoundary(    nElements*polynomialOrder,
                                                length,
                                                boundaryFlux) )

# Initial values
initialSolution     = np.zeros(model.size)

timeSeries          = solveNonlinearHeat1D( time,
                                            initialSolution,
                                            model,
                                            theta=finiteDifferenceImplicity,
                                            nonlinearSolver=nonlinearSolver,
                                            baseIncrement=baseIncrement,
                                            minIncrement=minIncrement,
                                            maxIncrement=maxIncrement,
                                            maxIncrements=maxIncrements,
                                            maxCorrections=maxCorrections,
                                            tolerance=tolerance,
                                            verbose=True,
                                            convergencePlot=convergencePlot )

# ---------------------------------------------------------
# Plot
animateTimeSeries(  time, 
                    samples, 
                    timeSeries, 
                    model,
                    speed=0.02,
                    ylim=( -0.1, 1.0 ) )

plt.tight_layout()
plt.show( block=True )