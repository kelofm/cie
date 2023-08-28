# --- Python imports ---
import numpy as np

# --- Visualization imports ---
from cie.glmesh import ParametricSurface, convertToSurfaceMesh, MeshNode, MeshApp3D

# --- Lighting imports ---
from cie.lighting import TimedSpotLight

# -----------------------------------------------------------
# Torus parameters
radius1     = 1.0
radius2     = 0.2

U           = np.linspace(0.0,2*np.pi,100)
V           = np.linspace(0.0,-2*np.pi,100)

torus   = ParametricSurface(    (   lambda u,v : np.cos(u)*(radius1+radius2*np.cos(v)),
                                    lambda u,v : np.sin(u)*(radius1+radius2*np.cos(v)),
                                    lambda u,v : radius2*np.sin(v) ),
                                domain=((min(U),max(U)),(min(V),max(V))) )
# -----------------------------------------------------------
# Pack geometry in a mesh node
geometry        = convertToSurfaceMesh( torus.evaluate(U,V) )
root            = MeshNode( vertices=geometry["vertices"], faces=geometry["faces"] )

# Create GL window
meshApp     = MeshApp3D(    root,
                            light=TimedSpotLight )

meshApp.view.camera.center          = np.asarray( [0.0,0.0,0.0], dtype=np.float32 )
meshApp.mesh.light._colorFunctor    = lambda t : (1.0,1.0,1.0)
meshApp.mesh.light._posFunctor      = lambda t : (0.0,0.0,10*(radius1+radius2)*np.cos(t))

meshApp.run()