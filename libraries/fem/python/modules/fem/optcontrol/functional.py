# --- Python Imports ---
import numpy as np

# --- Internal Imports ---
from cie.fem.numeric import Integrator

# ---------------------------------------------------------
def squaredSolutionErrorFunctional( solution, referenceSolution, model, integrationOrder=None ):
    '''
    Integrates the squared difference of the specified solution and a reference solution over the space domain.
    INTEGRATE{ (u_ref - u)^2 }dx

    Arguments:
        solution            : solution to compare the reference solution to
        referenceSolution   : desired solution
        model               : FEModel
        [integrationOrder]  : optional, default will use the double the first elements's integration order in the model
    '''
    # Initialize
    functionalValue = 0.0

    integrate = None
    if integrationOrder is None:
        integrate = Integrator( 2 * model.elements[0].integrator.polynomialOrder + 1 )
    else:
        integrate = Integrator( integrationOrder )
    
    for element in model.elements:
        localFunctional = lambda x: element( referenceSolution[element.DoFs] - solution[element.DoFs], x )**2
        functionalValue += integrate( localFunctional, element.domain )

    return functionalValue