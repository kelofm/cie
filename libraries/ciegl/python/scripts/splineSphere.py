# --- Splinekernel imports ---
from cie.pygeo.pysplinekernel import SurfaceKernel

# --- GL Visualization imports ---
from cie.lighting import SimpleLight
from cie.glmesh import convertToSurfaceMesh, MeshApp3D, MeshNode

# --- VisPy imports ---
from vispy.scene.visuals import create_visual_node
from vispy import scene

# --- Python imports ---
import numpy as np

# -----------------------------------------------------------

# Settings - sphere
radius                  = 1.0
offset                  = [ 0.0, 0.0, 0.0 ]
# Settings - spline
nInterpolationPoints    = [15,15]
polynomialOrders        = [3,3]
nSamples                = [10,10]

# -----------------------------------------------------------
# Define sphere
class Sphere:
        def __init__(self,radius=1.0,offset=[0.0,0.0,0.0]):
                self._radius    = radius
                if len(offset) is 2:
                        offset.append(0.0)
                if len(offset) is 3:
                        self._offset    = offset
                else:
                        raise ValueError('Invalid offset!')


        def surfacePoint(self, u, v):
                return [
                        self._radius*np.cos(v)*np.cos(u) - self._offset[0],
                        self._radius*np.cos(v)*np.sin(u) - self._offset[1],
                        self._radius*np.sin(v) - self._offset[2]
                ]


# Get points on the sphere (interpolation points)
S                       = Sphere(offset=[0.0,0.0,0.0])
interpolationPoints     = [[],[],[]]
for i,v in enumerate( np.linspace(np.pi/2,-np.pi/2,num=nInterpolationPoints[0]+1) ):
        interpolationPoints[0].append([])
        interpolationPoints[1].append([])
        interpolationPoints[2].append([])
        for j,u in enumerate( np.linspace(0.0,2*np.pi,num=nInterpolationPoints[1]+1) ):
                point   = S.surfacePoint(u,v)
                interpolationPoints[0][i].append(point[0])
                interpolationPoints[1][i].append(point[1])
                interpolationPoints[2][i].append(point[2])

# -----------------------------------------------------------
# Create spline surface
surf = SurfaceKernel(
    interpolationPoints=np.array(interpolationPoints),
    polynomialOrders=polynomialOrders)

# -----------------------------------------------------------
# Load geometry in a mesh node
geometry        = convertToSurfaceMesh( surf.generatePoints(nSamples=nSamples) )
root            = MeshNode( vertices=geometry["vertices"], faces=geometry["faces"] )

# Create GL window and load mesh
meshApp = MeshApp3D( root )
meshApp.view.camera.center = np.asarray( [0.0,0.0,0.0], dtype=np.float32 )
meshApp.mesh.light._pos = np.array( [0.0,0.0,-10.0], dtype=np.float32 )
meshApp.run()