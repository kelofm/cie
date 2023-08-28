# --- Python imports ---
import numpy as np
from timeit import default_timer

# --- GLMesh imports ---
from cie.glmesh import Sphere
from cie.glmesh import convertToSurfaceMesh
from cie.glmesh import MeshNode, MeshVisual, MeshApp3D

# --- GLTexture imports ---
from cie.gltexture import textureFolderPath, loadTexture
from cie.gltexture import uniformTexture, AggregatedTexture

# --- Vispy Imports ---
from vispy.app import Timer

# -----------------------------------------------------------
# Define geometries
sphere1     = Sphere(1.0,[0.0,0.0,0.0])
sphere2     = Sphere(1.0,[3.0,0.0,0.0])

# Sample
nSamples    = 50
sphere1     = convertToSurfaceMesh( sphere1.evaluate(sphere1.u(nSamples=nSamples), sphere1.v(nSamples=nSamples)) )
sphere2     = convertToSurfaceMesh( sphere2.evaluate(sphere2.u(nSamples=nSamples), sphere2.v(nSamples=nSamples)) )

# Load geometries in mesh nodes
root        = MeshNode( parent=None,
                        vertices=sphere1["vertices"],
                        faces=sphere1["faces"],
                        objectID=0,
                        materialID=0    )

root.addChild( MeshNode(    parent=root,
                            vertices=sphere2["vertices"],
                            faces=sphere2["faces"],
                            objectID=1,
                            materialID=0) )

# Get mesh sizes
meshSizes   = ( root.computeCompiledMeshSize(),
                root.child(0).computeCompiledMeshSize() )

# Create app
meshApp = MeshApp3D(root)
meshApp.view.camera.center  = [1.5,0.0,0.0]
meshApp.mesh.light._pos     = [1.5,0.0,0.0]

# Define translation function
t0 = default_timer()

def translateMesh(translation):
    for triangleIndex in range(len( root._vertices )):
        for componentIndex, displacementComponent in enumerate(translation):
            root._vertices[triangleIndex][componentIndex] = displacementComponent


def on_timer(event):
    t           = default_timer() - t0
    translation = [-np.sin(t), np.cos(t), 0.0]
    translateMesh(translation)
    print(root.get_vertices(indexed="faces")[0][0])
    meshApp.mesh._vertexBuffer.set_subdata( np.array(root.get_vertices(indexed="faces"), dtype=np.float32), offset=0 )
    meshApp.mesh.update()

timer = Timer(interval='auto', connect=on_timer, start='true')

# Run app
meshApp.run()