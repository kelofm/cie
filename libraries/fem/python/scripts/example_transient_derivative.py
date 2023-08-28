# --- Python Imports ---
import numpy as np
import scipy.sparse as sparse
from matplotlib import pyplot as plt

# --- Internal Imports ---
from cie.fem.discretization import IntegratedHierarchicBasisFunctions
from cie.fem.discretization import NonlinearHeatElement1D
from cie.fem.discretization import NonlinearTransientFEModel as Model
from cie.fem.discretization import DirichletBoundary
from cie.fem.postprocessing import ConvergencePlot

# ---------------------------------------------------------
# Geometry and material
length              = 1.0
#capacity            = lambda u: 1.0
#dCapacity           = lambda u: 0.0
capacity            = lambda u: 1.0 + 9.0 * np.exp( -(u-0.5)**2 / 0.005 )
dCapacity           = lambda u: 9.0 * np.exp( -(u-0.5)**2 / 0.005 ) * (2.0/0.005)*(0.5-u)
#conductivity        = lambda u: 1.0
#dConductivity       = lambda u: 0.0
conductivity        = lambda u: 1.0 + 9.0 * np.exp( -(u-0.5)**2 / 0.005 )
dConductivity       = lambda u: 9.0 * np.exp( -(u-0.5)**2 / 0.005 ) * (2.0/0.005)*(0.5-u)


# Load
load                = lambda t, x: 0.0
boundaryTemperature = 1.0

# Discretization
nElements           = 7
polynomialOrder     = 1

# Integration
integrationOrder    = 2 * (2*polynomialOrder + 1)
theta               = 0.5
dt                  = 0.1

# ---------------------------------------------------------
# General initialization
samples         = np.linspace( 0, length, num=100 )

# ---------------------------------------------------------
# Initialize FE model
model               = Model( nElements*polynomialOrder + 1, loadFunction=load )

# Create elements
basisFunctions      = IntegratedHierarchicBasisFunctions( polynomialOrder=polynomialOrder )
initialLoad         = lambda x: load(0,x)
model.elements      = [ NonlinearHeatElement1D( capacity,
                                                dCapacity,
                                                conductivity,
                                                dConductivity,
                                                (i*length/nElements, (i+1)*length/nElements),
                                                np.asarray( range(i*polynomialOrder, (i+1)*polynomialOrder+1) ),
                                                initialLoad,
                                                basisFunctions=basisFunctions,
                                                integrationOrder=integrationOrder   ) 
                        for i in range(nElements) ]

model.allocateZeros( )
derivative  = model.stiffness.copy()

# Boundary conditions (right BC is the load)
leftBCID    = model.addBoundaryCondition(   DirichletBoundary(  0, 
                                                                0.0,
                                                                lambda t: 0.0 ))

rightBCID   = model.addBoundaryCondition(   DirichletBoundary(  nElements*polynomialOrder,
                                                                length,
                                                                lambda t: boundaryTemperature) )

# ---------------------------------------------------------
u           = np.random.rand(model.size)
u[0]        = 0.0
u[-1]       = boundaryTemperature
du          = 1e-8
dr          = []

def residual( model, solution, loadFactor, massInverse=None ):
    if massInverse is None:
        massInverse = np.zeros(model.mass.shape, dtype=model.mass.dtype)
        sparse.linalg.inv(model.mass).todense( out=massInverse )
    temp = theta*massInverse.dot(model.stiffness.todense()) + 1.0/dt*np.eye(model.size)
    return np.ravel( temp.dot(u) - theta*massInverse.dot(model.load) )

def dResidual( model, solution, loadFactor, massInverse=None ):
        '''
        Compute the derivative of the residual with respect to the solution coefficients.
        See the readme for a better view: eq.(23)
        '''
        # Preallocate
        if massInverse is None:
            massInverse = np.zeros( model.mass.shape, dtype=model.mass.dtype )
            sparse.linalg.inv(model.mass).todense( out=massInverse )
        temp        = np.zeros(massInverse.shape)
        MKu         = massInverse.dot( model.stiffness.dot(solution) )
        Mq          = massInverse.dot(model.load)
        temperature = lambda x: model.sample(solution, x)

        # Integrate
        for e in model.elements:
            dCapacityCache = e.integrator.createCache(   lambda xi: e.capacityDerivative( temperature(e.toGlobalCoordinates(xi)) ),
                                                        e.basisFunctions.domain )
            for i in range(len(e.basisFunctions)):
                cache0 = dCapacityCache * e._basisCache[i]
                for l in range(len(e.basisFunctions)):
                    cache1 = cache0 * e._basisCache[l]
                    for j in range(len(e.basisFunctions)):
                        integral = e._jacobian * e.integrator.integrateCached(
                            lambda xi: 1.0,
                            cache1 * e._basisCache[j],
                            e.basisFunctions.domain
                            )
                        temp[e.DoFs[i], e.DoFs[l]] += integral * (loadFactor*Mq[e.DoFs[j]] - MKu[e.DoFs[j]])
        
        # Compute derivative
        return  theta * massInverse.dot(                            \
                            model.stiffness.todense()               \
                            + model.geometricStiffness.todense()    \
                            + temp                    )   \
                + np.eye(model.size)/dt

# Compute analytical derivative
model.resetMatrices()
model.integrate(    lambda x: model.sample(u, x),
                    u )
for boundary in model.boundaries:
    model.applyBoundaryCondition(boundary)

r               = residual(model, u, 1.0)
dra             = dResidual(model, u, 1.0)

# Compute finite difference derivative
for i in range(model.size):
    u[i]    += du

    # Recompute structural matrices
    model.resetMatrices()
    model.integrate(    lambda x: model.sample(u, x),
                        u )
    for boundary in model.boundaries:
        model.applyBoundaryCondition(boundary)

    diff = residual(model, u, 1.0) - r
    dr.append( diff/du )
    
    u[i]    -= du

# ---------------------------------------------------------
dr  = np.array(dr)
dr  = np.transpose( dr, (1,0) )
np.set_printoptions( precision=2, suppress=True )

print( dra - dr )
