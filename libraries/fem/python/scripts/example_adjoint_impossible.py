# --- Python Imports ---
import numpy as np
import scipy.sparse.linalg as linalg
from matplotlib import pyplot as plt

# --- Internal Imports ---
from cie.fem.numeric import solveLinearSystem
from cie.fem.discretization import IntegratedHierarchicBasisFunctions
from cie.fem.discretization import LinearHeatElement1D
from cie.fem.discretization import TransientFEModel
from cie.fem.discretization import DirichletBoundary, NeumannBoundary
from cie.fem.numeric import solveLinearHeat1D, solveAdjointLinearHeat1D
from cie.fem.optcontrol import squaredSolutionErrorFunctional
from cie.fem.postprocessing import animateTimeSeries, ConvergencePlot

'''
Usual 1D model but with a constant source on the farther half of the bar (impossible to correctly control with the given DoF)

The left boundary has a constant zero Dirichlet condition, while a constant zero Neumann condition is set on the right side.
The bar is subject to stationary source load, but only on the left half (top plot shows the source on the bar).

As always, the control is the flux on the right boundary, so following the applied load case is pretty much impossible, as
is shown by the results: the functional converges, but not to 0. Looking at the animated time history (middle plot), the control
overshoots on the right side, while not being able properly heat up the left half (as expected).
'''

# ---------------------------------------------------------
# SETTINGS
# ---------------------------------------------------------
# Reference
referenceControl            = lambda t: 0.0
initialControl              = lambda t: 0.0

# Geometry and material
length                      = 1.0
capacity                    = 1.0
conductivity                = 1.0

# Load
def load( t, x ):
    if x<0.5:
        return 1.0
    else:
        return 0.0

# Discretization
time                        = np.linspace(0.0, 1.0, 50)
nElements                   = 50
polynomialOrder             = 1

# Integration
integrationOrder            = 2*polynomialOrder + 1
finiteDifferenceImplicity   = 0.5

# Adjoint
numberOfAdjointIterations   = 25
regularization              = 10.0

# Postprocessing
numberOfSamples             = 100
convergencePlot             = ConvergencePlot()

# ---------------------------------------------------------
# REFERENCE SOLUTION
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
leftBCID    = model.addBoundaryCondition(   DirichletBoundary(  0,
                                                                0.0,
                                                                lambda t: 0.0   ) )

rightBCID   = model.addBoundaryCondition(   NeumannBoundary(    nElements*polynomialOrder,
                                                                length,
                                                                referenceControl) )

# Solve
model.updateTime( time[0] )
initialSolution     = np.zeros( model.size )
referenceTimeSeries = solveLinearHeat1D(    time,
                                            initialSolution,
                                            model,
                                            theta=finiteDifferenceImplicity )

# ---------------------------------------------------------
# OPTIMIZATION
# ---------------------------------------------------------
# Allocate
functionalValues    = np.zeros( numberOfAdjointIterations )
controls            = np.zeros( (numberOfAdjointIterations, len(time)) )

# Replace source
model._loadFunction = lambda t, x: 0.0

# Set initial control
u                   = np.asarray( [initialControl(t) for t in time] )

def controlFunction( t, control ):
    index = np.abs( time-t ).argmin()
    return control[index]

# Optimization loop
for i in range(numberOfAdjointIterations):

    # Set control
    model.boundaries[rightBCID].value   = lambda t: controlFunction( t, u )
    model.updateTime( time[0] )

    # Compute solution for the current control
    timeSeries  = solveLinearHeat1D(    time,
                                        initialSolution,
                                        model,
                                        theta=finiteDifferenceImplicity )

    initialAdjointSolution=None
    #########################################################################
    # Solve the stationary adjoint, and use the solution as the initial one
    initialAdjointSolution  = solveLinearSystem(    model.stiffness,
                                                    timeSeries[-1] - referenceTimeSeries[-1]    )
    #########################################################################

    # Compute adjoint solution
    adjointTimeSeries   = solveAdjointLinearHeat1D( time,
                                                    timeSeries - referenceTimeSeries,
                                                    model,
                                                    theta=finiteDifferenceImplicity,
                                                    initialAdjointSolution=initialAdjointSolution   )

    # Compute projections and update control
    projections         = -1.0/regularization * np.asarray([
                                adjointTimeSeries[tIndex][-1] for tIndex, t in enumerate(time)
                                ])
    #projections         = np.asarray([ np.min((0.04,np.max((0.0,projection)))) for projection in projections ])
    u                   += projections

    # Get functional integrated in space (for every time step)
    functionalValue = [ squaredSolutionErrorFunctional( solution,
                                                        referenceSolution,
                                                        model   )
                        for solution, referenceSolution in zip(timeSeries, referenceTimeSeries) ]

    # Integrate functional in time
    functionalValue = np.trapz( functionalValue,
                                x=time  )
    convergencePlot( functionalValue )

    # Save intermediate results
    controls[i]         = u.copy()
    functionalValues[i] = functionalValue

    # Print detauls
    message     = "Iteration\t\t:%i" % i
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
model.boundaries[rightBCID].value   = lambda t: controlFunction( t, u )
timeSeries = solveLinearHeat1D( time,
                                initialSolution,
                                model,
                                theta=finiteDifferenceImplicity )


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
axes[3].plot( u, ".-" )
axes[3].set_xlabel( "time" )
axes[3].set_ylabel( "Control" )

fig.canvas.draw()
plt.show( block=True )