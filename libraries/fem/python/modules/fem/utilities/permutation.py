# --- Python Imports ---
import itertools
import numpy as np

# ---------------------------------------------------------
def binaryPermutation( numberOfBits, dtype=np.int ):
    return np.asarray( [list(i) for i in itertools.product([0,1], repeat=numberOfBits)], dtype=dtype )


def oneHotArray( numberOfBits, dtype=np.int ):
    return np.asarray( [ [int(i==j) for i in range(numberOfBits)] for j in range(numberOfBits) ], dtype=dtype )


def gridPairs( components ):
    return np.asarray( [list(i) for i in itertools.product( components, repeat=2 )] )