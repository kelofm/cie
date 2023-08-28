# --- Python imports ---
import numpy as np

# -----------------------------------------------------

def convertToSurfaceMesh(grid):
    size0 = len(grid[0])
    size1 = len(grid[0][0])
    geometry = {'vertices':[], 'faces':[]}
    geometry['vertices'] = np.reshape(
                    np.transpose( 
                        np.array( grid ,dtype=np.float32  ), 
                        axes=(1,2,0)),
                    ( size0*size1 ,3))
    geometry['faces'] = np.zeros( ( 2*(size0-1)*(size1-1), 3 ), dtype=np.uint32 )
    k = 0
    for j in range(size0-1):
        for i in range(size1-1):
            jni = j*size1+i
            geometry['faces'][k] = [ 
                jni, 
                jni+size1,
                jni+size1+1
                ]
            geometry['faces'][k+1] = [ 
                jni, 
                jni+1+size1,
                jni+1
                ]
            k+=2
    return geometry