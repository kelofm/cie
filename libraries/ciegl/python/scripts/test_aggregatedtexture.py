# --- Python imports ---
import numpy as np

# --- Texture imports ---
from cie.gltexture import loadTexture, emptyTexture
from cie.gltexture import AggregatedTexture

# --- Vispy Imports ---
from vispy.gloo import Texture2D

# -----------------------------------------------------------
agg     = AggregatedTexture(shape=(10,10))
texture = np.asarray(   [   [ (0.0,0.0,0.0),(0.5,0.5,0.5)],
                            [ (0.5,0.5,0.5),(0.0,0.0,0.0)]   ],
                        dtype=np.float32 )

agg.registerTexture( emptyTexture( (5,5) ) ,2)
agg.registerTexture(texture,0)
agg.registerTexture(texture,1)

vertexShader, fragmentShader = agg.generateShaderCode()

with open('vertexShader.txt','w') as file:
    file.write(vertexShader)

with open('fragmentShader.txt','w') as file:
    file.write(fragmentShader)