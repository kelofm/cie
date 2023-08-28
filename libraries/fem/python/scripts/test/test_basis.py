# --- Python Imports ---
import numpy as np
import matplotlib.pyplot as plt

# --- PyFEM Imports ---
from cie.fem.utilities import *
from cie.fem.discretization import *
from cie.fem.numeric import *
from cie.fem.discretization.basisfunctions import BasisFunctions

# ---------------------------------------------------------
testPolynomial  = lambda t: 1.0 + 2.0*t + 3.0*t*t + 4.0*t*t*t
testIntegral    = lambda t: t*(1.0 + t*(1.0 + t*(1.0 + t)))
integrator      = Integrator(3)

domain = (-1.0,1.0)
testValue   = testIntegral(domain[1]) - testIntegral(domain[0])
value = integrator( testPolynomial, domain )
if np.abs(value-testValue) < 1e-15:
    print("PASS")
else:
    print("FAIL " + str(value))

value = integrator( testPolynomial, (0.0, 1.0) )
if np.abs(value-4.0) < 1e-15:
    print("PASS")
else:
    print("FAIL " + str(value))

value = integrator( testPolynomial, (0.0, 5.0) )
if np.abs(value-780.0) < 1e-15:
    print("PASS")
else:
    print("FAIL " + str(value))

polynomialOrder = 3
basisFunctions  = IntegratedHierarchicBasisFunctions( polynomialOrder=polynomialOrder )
#basisFunctions  = HierarchicBasisFunctions( polynomialOrder=polynomialOrder )
#basisFunctions  = BasisFunctions( functions=[np.polynomial.legendre.Legendre(coefficients) for coefficients in oneHotArray(polynomialOrder+1)] )
t               = np.linspace( 2*basisFunctions.domain[0], 2*basisFunctions.domain[1], num=100 )

print( "Nodal functions:\t" + str(basisFunctions.nodalFunctions) )

for basisID in range( polynomialOrder+1 ):
    plt.plot( t, basisFunctions( basisID, t ) )
    integral = np.polynomial.legendre.legint( basisFunctions._functions[basisID].coef )
    integral = np.polynomial.legendre.Legendre( integral )
    print( integral(1.0)-integral(-1.0) )

#derivatives     = basisFunctions.derivatives()
#for basisID in range( polynomialOrder+1 ):
#    plt.plot( t, derivatives( basisID, t ) )


plt.show()