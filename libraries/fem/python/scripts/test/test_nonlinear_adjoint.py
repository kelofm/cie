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

from cie.fem.optcontrol.adjointmodel import AdjointHeatElement1D, AdjointModel
from cie.fem.numeric import solveAdjointNonlinearHeat1D
from cie.fem.numeric import solveLinearSystem
from cie.fem.optcontrol import squaredSolutionErrorFunctional

# ---------------------------------------------------------
# Geometry
length                      = 1.0

# Material
a                           = 1.0
b                           = 9.0
c                           = 0.5
d                           = 0.005

exp                         = lambda u: np.exp( -(u-c)**2 / d )

capacity                    = lambda u: a + b * exp(u)
dCapacity                   = lambda u: b * exp(u) * (2.0/d)*(c-u)

conductivity                = capacity
dConductivity               = dCapacity
ddConductivity              = lambda u: 2.0*b/d * exp(u) * ( 2.0/d*(c-u)**2 - 1.0 )

# Load
load                        = lambda t, x: 0.0
penaltyValue                = 1e10

#referenceControl            = lambda t: 2.0*(1.0 - np.exp(-100.0*t)) - t/2.0
referenceControl            = lambda t: t
initialControl              = lambda t: 0.0

# Adjoint settings
numberOfAdjointIterations   = 10
regularization              = 1e1

# Discretization
time                        = np.linspace(0.0, 1.0, num=25)
nElements                   = 10
polynomialOrder             = 1

# Integration
integrationOrder            = 1 * (2*polynomialOrder + 1)
finiteDifferenceImplicity   = 0.5

# Iteration
baseIncrement       = 1.0
minIncrement        = 0.01
maxIncrement        = 1.0
maxIncrements       = 35
maxCorrections      = 15
tolerance           = 1e-10

# Postprocessing and visualization
numberOfSamples     = 100

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
leftBCID    = model.addBoundaryCondition( DirichletBoundary(    0, 
                                                                0.0,
                                                                lambda t: 0.0   ) )

rightBCID   = model.addBoundaryCondition( NeumannBoundary(  nElements*polynomialOrder,
                                                            length,
                                                            referenceControl) )

# Initial values
initialSolution     = np.zeros(model.size)


# ---------------------------------------------------------
# REFERENCE SOLUTION
# ---------------------------------------------------------
referenceTimeSeries = solveNonlinearHeat1D( time,
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
                                            verbose=False,
                                            convergencePlot=convergencePlot )

# ---------------------------------------------------------
# APPROXIMATE GRADIENT
# ---------------------------------------------------------
dControl            = 1e-8
referenceControls   = np.array( [ model.boundaries[rightBCID].value(t) for t in time ] )
initialControls     = np.array( [ initialControl(t) for t in time ] )

def makeControl( index ):
    def control( t ):
        return initialControls[index] + dControl
    return control

gradient    = []

for i in range(len(time)):
    model.boundaries[rightBCID].value = makeControl(i)
    timeSeries = solveNonlinearHeat1D(  time,
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
                                        verbose=False,
                                        convergencePlot=convergencePlot )
    # Get functional integrated in space (for every time step)
    functionalValue = [ squaredSolutionErrorFunctional( solution,
                                                        referenceSolution,
                                                        model   )
                        for solution, referenceSolution in zip(timeSeries, referenceTimeSeries) ]

    # Integrate functional in time
    functionalValue = np.trapz( functionalValue,
                                x=time  )
    gradient.append( functionalValue/dControl )

gradient    = np.asarray(gradient)
gradient    = gradient / np.linalg.norm(gradient)

# ---------------------------------------------------------
# ADJOINT GRADIENT
# ---------------------------------------------------------
adjointModel            = AdjointModel( model,
                                        time,
                                        referenceTimeSeries,
                                        timeSeries,
                                        ddConductivity )
adjointModel.allocateZeros()
adjointModel.addBoundaryCondition(  model.boundaries[leftBCID] )
adjointModel.addBoundaryCondition(  model.boundaries[rightBCID] )

# Reset adjoint BC (control)
adjointModel.boundaries[rightBCID].value = initialControl

#########################################################################
# Solve the stationary adjoint, and use the solution as the initial one
adjointModel.updateTime( len(time)-1 )
initialAdjointSolution  = solveLinearSystem(    adjointModel.stiffness + adjointModel.nonsymmetricStiffness,
                                                timeSeries[-1] - referenceTimeSeries[-1]    )
#########################################################################

adjointTimeSeries       = solveAdjointNonlinearHeat1D(  adjointModel,
                                                        theta=finiteDifferenceImplicity,
                                                        initialAdjointSolution=None )

# Compute projections
projections         = -1.0/regularization * np.asarray([ 
                            adjointTimeSeries[tIndex][-1] for tIndex, t in enumerate(time)
                            ])
projections = projections / np.linalg.norm(projections)
# ---------------------------------------------------------
# ---------------------------------------------------------
print( gradient )
print( projections )
print( gradient-projections )