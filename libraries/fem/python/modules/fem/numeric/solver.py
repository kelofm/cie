# --- Python Imports ---
import numpy as np
import scipy
import scipy.sparse.linalg as linalg

# ---------------------------------------------------------
def solveLinearSystem( A, b, *args, x0=None, sparse=True, **kwargs ):
    '''
    Wrapper function to simplify changing linear solvers in the future
    '''
    # Initialize
    if x0 is None:
        x0 = np.zeros( b.shape )

    info = 0

    if sparse:
        # Iterative solver
        #sol, info   = scipy.sparse.linalg.gmres(   A,
        #                                           b,
        #                                           x0=x0,
        #                                           atol=1e-12,
        #                                           maxiter=np.max((5*len(b),100)) )

        # Direct solver
        sol     = scipy.sparse.linalg.spsolve(  A,
                                                b )

    else:
        # Iterative solver
        #sol, info   = scipy.linalg.gmres(  A,
        #                                   b,
        #                                   x0=x0,
        #                                   atol=1e-12,
        #                                   maxiter=np.max((5*len(b),100)) )

        # Direct solver
        sol     = scipy.linalg.solve(   A,
                                        b )

    if info < 0:
        raise RuntimeError( "Linear solver error!" )
    elif info>0:
        print( "Solver failed to converge after " + str(info) + " steps" )

    return sol