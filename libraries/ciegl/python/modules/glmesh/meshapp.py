# --- Python imports ---
import numpy as np

# --- Internal imports ---
from . import MeshVisual
from . import defaultVertexShader, defaultFragmentShader

# --- Texture imports ---
from cie.gltexture import uniformTexture, AggregatedTexture

# --- Lighting imports ---
from cie.lighting import SimpleLight

# --- Vispy imports ---
from vispy import scene, app

# -----------------------------------------------------
class MeshApp3D:
    def __init__(self,
                        meshRoot,
                        aggregatedTexture=None,
                        camera=None,
                        bgcolor=(0.2,0.2,0.2),
                        border_color=(0.0,0.0,0.0),
                        border_width=1,
                        **kwargs ):
        # Default texture
        if aggregatedTexture is None:
            aggregatedTexture   = AggregatedTexture(shape=(2,2))
            aggregatedTexture.registerTexture(  uniformTexture(size=(2,2)),
                                                0,
                                                ambientMaterialConstant=0.7,
                                                specularMaterialConstant=0.8    )

        # Vispy init
        self._meshNode          = scene.visuals.create_visual_node(MeshVisual)
        self._canvas            = scene.SceneCanvas(    keys='interactive',
                                                        size=(1024,768) )
        self.view               = self._canvas.central_widget.add_view(     bgcolor=bgcolor,
                                                                            border_color=border_color,
                                                                            border_width=border_width )

        # Default arguments
        if camera is None:
            self.view.camera   = scene.cameras.ArcballCamera(  fov=0,
                                                                center=(0.0,0.0,0.0),
                                                                parent=self.view )
        else:
            self.view.camera   = camera( parent=self.view )

        # Create visual
        self.mesh               = self._meshNode(   meshRoot,
                                                    aggregatedTexture,
                                                    camera=self.view.camera,
                                                    **kwargs )

        # Adjust lighting
        vertices                    = self.mesh.root.getCompiledMesh()['vertices']
        vertices                    = np.reshape( vertices, (len(vertices)*3,3) )

        center                      = ( ( np.max(vertices[:,0]) - np.min(vertices[:,0]) )/2.0,
                                        ( np.max(vertices[:,1]) - np.min(vertices[:,1]) )/2.0,
                                        ( np.max(vertices[:,2]) - np.min(vertices[:,2]) )/2.0     )

        self.mesh.light._pos        = center
        self.view.camera.center     = center
        self.mesh.light._color      = np.array( (1.0,1.0,1.0), dtype=np.float32 )

        self.view.add(self.mesh)



    def run(self):
        self._canvas.show()
        app.run()