# --- Python Imports ---
import numpy as np

# ---------------------------------------------------------

def isArray( var ):
    if type(var) in (list, tuple, np.ndarray):
        return True
    else:
        return False


def isNumpyArray( var ):
    if type(var) is np.ndarray:
        return True
    else:
        return False