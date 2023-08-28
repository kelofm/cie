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
numberOfAdjointIterations   = 150
regularization              = 5e0

# Discretization
time                        = np.linspace(0.0, 1.0, num=20)
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
# ADJOINT LOOP
# ---------------------------------------------------------
referenceControl    = np.array( [ model.boundaries[rightBCID].value(t) for t in time ] )
functionalValues    = np.zeros( numberOfAdjointIterations )
controls            = np.zeros( (numberOfAdjointIterations, len(time)) )

# Set initial control
control             = np.array( [initialControl(t) for t in time] )

def controlFunctor( t, control ):
    index = np.abs( time-t ).argmin()
    return control[index]

# Optimization loop
for i in range(numberOfAdjointIterations):

    # Set control
    model.boundaries[rightBCID].value   = lambda t: controlFunctor(t,control)
    model.updateTime( time[0], initialSolution )

    # Compute solution for the current control
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

    adjointModel            = AdjointModel( model,
                                            time,
                                            referenceTimeSeries,
                                            timeSeries,
                                            ddConductivity )
    adjointModel.allocateZeros()
    adjointModel.addBoundaryCondition(  model.boundaries[leftBCID] )
    adjointModel.addBoundaryCondition(  model.boundaries[rightBCID] )

    initialAdjointSolution = None
    #########################################################################
    # Solve the stationary adjoint, and use the solution as the initial one
    adjointModel.updateTime( len(time)-1 )
    initialAdjointSolution  = solveLinearSystem(    model.stiffness,
                                                    timeSeries[-1] - referenceTimeSeries[-1]    )
    #########################################################################

    adjointTimeSeries       = solveAdjointNonlinearHeat1D(  adjointModel,
                                                            theta=finiteDifferenceImplicity,
                                                            initialAdjointSolution=initialAdjointSolution )

    # Compute projections and update control
    projections         = -1.0/regularization * np.asarray([ 
                                adjointTimeSeries[tIndex][-1] for tIndex, t in enumerate(time)
                                ])
    control             += projections

    # Get functional integrated in space (for every time step)
    functionalValue = [ squaredSolutionErrorFunctional( solution,
                                                        referenceSolution,
                                                        model   )
                        for solution, referenceSolution in zip(timeSeries, referenceTimeSeries) ]

    # Integrate functional in time
    functionalValue = np.trapz( functionalValue,
                                x=time  )

    # Save intermediate results
    controls[i]         = control.copy()
    functionalValues[i] = functionalValue

    # Print detauls
    message     = "Iteration\t\t:%i" % i
    message     += "\nNorm of control error\t:%.3E" % np.trapz( (referenceControl - control)**2, x=time )
    message     += "\nFunctional value\t:%.3E" % functionalValue
    message     += "\n"
    print( message )


# ---------------------------------------------------------
# POSTPROCESSING
# ---------------------------------------------------------
# Find best control
uIndex  = functionalValues.argmin()
u       = controls[uIndex]

# Compute solution with best control
model.boundaries[rightBCID].value   = lambda t: controlFunctor( t, u )
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
                                    verbose=False )


# Output and graphics
print( "\nFinal control:\t" + str(u) + "\nin iteration " + str(uIndex) )

fig     = plt.figure( )
axes    = ( fig.add_subplot( 3,1,1 ),
            fig.add_subplot( 3,1,2 ),
            fig.add_subplot( 3,2,5 ),
            fig.add_subplot( 3,2,6 ))

# Load
samples             = np.linspace( 0, length, num=numberOfSamples )
axes[0].plot( samples, [ load(0,x) for x in samples ], "r--" )
axes[0].set_xlabel( "x" )
axes[0].set_ylabel( "load" )

# Animation
animateTimeSeries(  time,
                    samples,
                    np.asarray((referenceTimeSeries, timeSeries)),
                    model,
                    speed=0.05,
                    ylim=(-0.1,0.5),
                    figure=fig,
                    axis=axes[1] )
axes[1].set_xlabel( "x" )
axes[1].set_ylabel( "y" )
axes[1].legend( ["Reference solution", "Final solution", "Load of reference"] )

# Functional
axes[2].plot( functionalValues, ".-" )
axes[2].set_yscale( "log" )
axes[2].set_xlabel( "# Iteration" )
axes[2].set_ylabel( "Functional value" )
axes[2].grid( b=True, axis="y" )

# Final control
axes[3].plot( referenceControl, ".-" )
axes[3].plot( u, ".-" )
axes[3].set_xlabel( "time" )
axes[3].set_ylabel( "Control" )
axes[3].legend( [ "ref", "opt" ] )

fig.canvas.draw()
plt.show( block=True )


# Write to file
with open("controls.csv", "w") as file:
    for i, control in enumerate(controls):
        file.write( str(i)+"," )
        for value in control:
            file.write( str(value)+"," )
        file.write("\n")