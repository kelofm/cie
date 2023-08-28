# --- Python Imports ---
import numpy as np
import scipy.sparse.linalg as linalg
from matplotlib import pyplot as plt

# --- Internal Imports ---
from cie.fem.discretization import IntegratedHierarchicBasisFunctions
from cie.fem.discretization import LinearHeatElement1D
from cie.fem.discretization import TransientFEModel
from cie.fem.discretization import DirichletBoundary, NeumannBoundary
from cie.fem.numeric import solveLinearHeat1D
from cie.fem.optcontrol import squaredSolutionErrorFunctional

# ---------------------------------------------------------
# Geometry and material
length                      = 1.0
capacity                    = 1.0
conductivity                = 1.0

# Load
load                        = lambda t, x: 0.0

# Discretization
time                        = np.linspace(0.0, 1.0, 25)
nElements                   = 25
polynomialOrder             = 3

# Integration
integrationOrder            = 2*polynomialOrder + 1
finiteDifferenceImplicity   = 0.5

# ---------------------------------------------------------
# Initialize FE model
model               = TransientFEModel( nElements*polynomialOrder + 1, loadFunction=load )

# Create elements
basisFunctions      = IntegratedHierarchicBasisFunctions( polynomialOrder=polynomialOrder )
initialLoad         = lambda x: load( time[0], x )
model.elements      = [ LinearHeatElement1D(    capacity,
                                                conductivity,
                                                (i*length/nElements, (i+1)*length/nElements),
                                                np.asarray( range(i*polynomialOrder, (i+1)*polynomialOrder+1) ),
                                                initialLoad,
                                                basisFunctions=basisFunctions,
                                                integrationOrder=integrationOrder   )
                        for i in range(nElements) ]

# Integrate
model.allocateZeros( )
model.integrate( )

# Boundary conditions
penaltyValue    = 1e3
leftBCID    = model.addBoundaryCondition(   DirichletBoundary(  0,
                                                                0.0,
                                                                lambda t: 0.0,
                                                                penaltyValue=penaltyValue   ) )

rightBCID   = model.addBoundaryCondition(   NeumannBoundary(    nElements*polynomialOrder,
                                                                length,
                                                                lambda t: 1.0) )

# Solve
initialSolution     = np.zeros( model.size )
referenceTimeSeries = solveLinearHeat1D(    time,
                                            initialSolution,
                                            model,
                                            theta=finiteDifferenceImplicity )

# ---------------------------------------------------------
# Modify boundary conditions
for neumannBCValue in 1.0 - np.linspace(0.0,1.0,11):
    model.boundaries[rightBCID].value = lambda t: neumannBCValue

    # Solve new system
    timeSeries      = solveLinearHeat1D(    time,
                                            initialSolution,
                                            model,
                                            theta=finiteDifferenceImplicity )

    # Get functional integrated in space
    functionalValue = [ squaredSolutionErrorFunctional( solution,
                                                        referenceSolution,
                                                        model   )
                        for solution, referenceSolution in zip(timeSeries, referenceTimeSeries) ]

    # Integrate in time
    functionalValue = np.trapz( functionalValue,
                                x=time  )

    print( ("Boundary load = %.1f" % neumannBCValue) + ",\tFunctional = %.5f" % functionalValue )