# --- Python Imports ---
import numpy as np
import scipy.sparse.linalg as linalg
from matplotlib import pyplot as plt

# --- Internal Imports ---
from cie.fem.numeric import solveLinearSystem
from cie.fem.discretization import IntegratedHierarchicBasisFunctions as Basis
from cie.fem.discretization import NonlinearHeatElement1D as Element
from cie.fem.discretization import NonlinearTransientFEModel as Model
from cie.fem.discretization import DirichletBoundary, NeumannBoundary
from cie.fem.numeric import solveNonlinearHeat1D
from cie.fem.optcontrol import squaredSolutionErrorFunctional
from cie.fem.postprocessing import animateTimeSeries, ConvergencePlot

# ---------------------------------------------------------
# SETTINGS
# ---------------------------------------------------------
# Reference
referenceControl            = lambda t: 1.0
initialControl              = lambda t: 0.0

# Geometry and material
length                      = 1.0
capacity                    = lambda u: 1.0 + 9.0 * np.exp( -(u-0.5)**2 / 0.005 )
dCapacity                   = lambda u: 9.0 * np.exp( -(u-0.5)**2 / 0.005 ) * (2.0/0.005)*(0.5-u)
conductivity                = lambda u: 1.0 + 9.0 * np.exp( -(u-0.5)**2 / 0.005 )
dConductivity               = lambda u: 9.0 * np.exp( -(u-0.5)**2 / 0.005 ) * (2.0/0.005)*(0.5-u)

# Load
load                        = lambda t, x: 0.0

# Discretization
time                        = np.linspace(0.0, 1.0, 15)
nElements                   = 50
polynomialOrder             = 1
penaltyValue                = 1e10

# Integration
integrationOrder            = 2*polynomialOrder + 1
finiteDifferenceImplicity   = 0.5

# Adjoint
numberOfAdjointIterations   = 50
regularization              = 8.0

# Postprocessing
numberOfSamples             = 100
convergencePlot             = ConvergencePlot()

# ---------------------------------------------------------
# REFERENCE SOLUTION
# ---------------------------------------------------------
# Initialize FE model
model               = Model( nElements*polynomialOrder + 1, loadFunction=load )

# Create elements
basisFunctions      = Basis( polynomialOrder=polynomialOrder )
initialLoad         = lambda x: load( time[0], x )
model.elements      = [ Element(capacity,
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
                                                                lambda t: 0.0,
                                                                penaltyValue=penaltyValue   ) )

rightBCID   = model.addBoundaryCondition(   NeumannBoundary(    nElements*polynomialOrder,
                                                                length,
                                                                referenceControl) )

# Solve
model.updateTime( time[0] )
initialSolution     = np.zeros( model.size )
referenceTimeSeries = solveNonlinearHeat1D(    time,
                                            initialSolution,
                                            model,
                                            theta=finiteDifferenceImplicity )

# ---------------------------------------------------------
# OPTIMIZATION
# ---------------------------------------------------------
# Allocate
functionalValues    = np.zeros( numberOfAdjointIterations )
controls            = np.zeros( (numberOfAdjointIterations, len(time)) )
referenceControl    = np.asarray( [ model.boundaries[rightBCID].value(t) for t in time ] )

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
    timeSeries  = solveNonlinearHeat1D(    time,
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
    message     += "\nNorm of control error\t:%.3E" % np.trapz( (referenceControl - u)**2, x=time )
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
axes    = ( fig.add_subplot( 2,1,1 ),
            fig.add_subplot( 2,2,3 ),
            fig.add_subplot( 2,2,4 ))

# Animation
samples             = np.linspace( 0, length, num=numberOfSamples )
animateTimeSeries(  time,
                    samples,
                    np.asarray((referenceTimeSeries, timeSeries)),
                    model,
                    speed=0.05,
                    ylim=(-0.1,1.0),
                    figure=fig,
                    axis=axes[0] )
axes[0].set_xlabel( "x" )
axes[0].set_ylabel( "y" )
axes[0].legend( ["Reference solution", "Final solution"] )

# Functional
axes[1].plot( functionalValues, ".-" )
axes[1].set_yscale( "log" )
axes[1].set_xlabel( "# Iteration" )
axes[1].set_ylabel( "Functional value" )
axes[1].grid( b=True, axis="y" )

# Final error in control
axes[2].plot( referenceControl - u, ".-" )
axes[2].set_xlabel( "time" )
axes[2].set_ylabel( "Control error" )

fig.canvas.draw()
plt.show( block=True )