# --- Python Imports ---
import numpy as np

# --- PyFEM Imports ---
from cie.fem.numeric import *

# ---------------------------------------------------------
def test( function, primitiveFunction, integrator, domain ):
    # Test integration
    value       = integrator( function, domain )    
    testValue   = primitiveFunction(domain[1]) - primitiveFunction(domain[0])
    error       = np.abs( value-testValue )
    
    if error < 1e-15:
        print("PASS")
    else:
        print("FAIL " + str(error))

    # Test cached integration
    cache       = integrator.createCache( function, domain )
    cachedValue = integrator.integrateCached( lambda x: 1.0, cache, domain )
    error       = np.abs( value-cachedValue )

    if error > 1e-16:
        print( "Cached integration failed! " + str(error) )

polynomialOrder = 3

coefficients    = np.random.rand( polynomialOrder+1 )
#coefficients    = np.arange( 1, polynomialOrder+2 )

def testPolynomial( t, order=0, coefficients=coefficients ):
    if order<polynomialOrder:
        return coefficients[order] + t * testPolynomial( t, order=order+1 )
    else:
        return coefficients[order]

def temp( t, order=1, coefficients=coefficients ):
    if order<polynomialOrder+1:
        return coefficients[order-1]/float(order) + t * temp(t, order=order+1)
    else:
        return coefficients[order-1]/float(order)

testIntegral = lambda t: t*temp(t)

# ---------------------------------------------------------
integrator      = Integrator( polynomialOrder )
testDomains     = ( (-1.0, 1.0),
                    (0.0, 1.0),
                    (1.0, 5.0),
                    (0.0, 100.0),
                    (0.1, 0.11),
                    (np.pi, np.pi**5)   )

for domain in testDomains:
    test(   testPolynomial,
            testIntegral,
            integrator,
            domain )