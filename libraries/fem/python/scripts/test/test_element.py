# --- Python Imports ---
import numpy as np
from matplotlib import pyplot as plt

# --- Internal Imports ---
from cie.fem.discretization import IntegratedHierarchicBasisFunctions
from cie.fem.discretization import LinearHeatElement1D

# ---------------------------------------------------------

tolerance       = 1e-4

length          = 3.14159265

a               = 1.0
b               = 9.0
c               = 0.5
d               = 0.005
exp             = lambda u: np.exp( -(u-c)**2 / d )

capacity        = lambda u: a + b * exp(u)
dCapacity       = lambda u: b * exp(u) * (2.0/d)*(c-u)
conductivity    = capacity
dConductivity   = dCapacity

polynomialOrder = 3
basisFunctions  = IntegratedHierarchicBasisFunctions( polynomialOrder=polynomialOrder )
element         = LinearHeatElement1D(  capacity,
                                        conductivity,
                                        [0.0,length],
                                        [i for i in range(polynomialOrder+1)],
                                        lambda x: 0.0,
                                        basisFunctions=basisFunctions,
                                        integrationOrder=polynomialOrder )

# Test basis function derivatives
dx                      = 1e-12
samples                 = np.linspace(0.0,length-dx,num=100)
for function, derivative in zip(element.basisFunctions, element.basisDerivatives):
    approximateDerivative   = lambda x: ( function(x+dx) - function(x) ) / dx
    squaredDifference       = lambda x: (approximateDerivative(x)-derivative(x))**2
    err                     = np.trapz( [squaredDifference(x) for x in samples],
                                        x=samples )
    print( "Basis derivative error: %.3f" % err )
    assert( err < tolerance )

# Test derivative
u                       = np.random.rand(polynomialOrder+1)
approximateDerivative   = lambda x: ( element(u,x+dx) - element(u,x) ) / dx
derivative              = lambda x: element.derivative(u,x)
squaredDifference       = lambda x: (approximateDerivative(x)-derivative(x))**2

err                     = np.trapz( [squaredDifference(x) for x in samples],
                                    x=samples )

print( "Error: %.5f" % err )

#if err > 1e-3:
#    plt.plot( samples, element(u,samples), "k.-" )
#    plt.plot( samples, [approximateDerivative(x) for x in samples], "b" )
#    plt.plot( samples, [derivative(x) for x in samples], "r" )
#    plt.show()

assert( err < tolerance )