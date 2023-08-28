# --- Python imports ---
from sys import getsizeof as sizeof
import numpy as np

# --- GLMesh imports ---
from cie.glmesh import MeshNode, MeshApp3D

# --- GLTexture imports ---
from cie.gltexture import textureFolderPath, loadTexture
from cie.gltexture import AggregatedTexture

# -----------------------------------------------------------
v1      = np.asarray(   [   [0.0,0.0,0.0],
                            [1.0,0.0,0.0],
                            [1.0,1.0,0.0],
                            [0.0,1.0,0.0]],
                        dtype=np.float32)
f1      = np.asarray(   [   [0,1,2],
                            [0,2,3]],
                            dtype=np.uint32)
tc1     = np.asarray(   [   [0.0,0.0],
                            [1.0,0.0],
                            [1.0,1.0],
                            [0.0,1.0]],
                            dtype=np.float32)

# -----------------------------------------------------------
textureXY   = loadTexture(textureFolderPath + '/cage.png')
textureDim  = ( 1.5*textureXY.shape[0], 1.5*textureXY.shape[1] )

agg         = AggregatedTexture( shape=textureDim )
agg.registerTexture( textureXY, 0, ambientMaterialConstant=0.5, diffuseMaterialConstant=0.5, specularMaterialConstant=1.0 )

# -----------------------------------------------------------
root        = MeshNode( vertices=v1,
                        faces=f1,
                        textureCoordinates=tc1,
                        objectID=0,
                        materialID=0 )

# -----------------------------------------------------------
mApp    = MeshApp3D(    root,
                        agg )

mApp.mesh.light._pos    = np.array((0.5,0.5,0.5),dtype=np.float32)

floatSize       = 4
pointSize       = 3*floatSize
triangleSize    = 3*pointSize

print(mApp.mesh._vertexBuffer.nbytes)
newPos          = np.array( ((0.0,-1.0,0.0),(1.0,-1.0,0.0),(1.0,0.0,0.0)), dtype=np.float32 )
mApp.mesh._vertexBuffer.set_subdata(newPos,offset=3)

mApp.run()