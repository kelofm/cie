from vispy import scene, app

import numpy as np
from cie.glmesh import TriangleMeshVisual

from cie.gltexture import defaultFragmentShaderWithTexture, defaultVertexShaderWithTexture

vertices    = np.array( [   (0.0, 0.0, 0.0),
                            (1.0, 0.0, 0.0),
                            (1.0, 0.0, 1.0),
                            (0.0, 0.0, 1.0) ],
                        dtype=np.float32 )

faces       = np.array( [   (0,3,2),
                            (0,2,1)],
                        dtype=np.uint32)

texture     = np.array( [   ( (1.0,0.0,0.0), (0.0,1.0,0.0) ),
                            ( (0.0,0.0,1.0), (1.0,1.0,1.0) ) ],
                        dtype=np.float32 )


textureCoordinates  = np.array( [   (0.0, 0.0),
                                    (1.0, 0.0),
                                    (1.0, 1.0),
                                    (0.0, 1.0) ],
                                dtype=np.float32)


# Add mesh as VisPy node
TriangleMesh = scene.visuals.create_visual_node(TriangleMeshVisual)


# Create a new OpenGL window
plt = scene.SceneCanvas(        keys='interactive',
                                size=(1024, 768) )

view = plt.central_widget.add_view(     bgcolor=(0.2,0.2,0.2),
                                        border_color=(0.0,0.0,0.0),
                                        border_width=1 )

# Add and configure camera
view.camera = scene.cameras.ArcballCamera(fov=0)
view.camera.center  = (0.5,0.0,0.5)

mesh    = TriangleMesh(             vertices,
                                    faces,
                                    texture=texture,
                                    textureCoordinates=textureCoordinates,
                                    vertexShader=defaultVertexShaderWithTexture,
                                    fragmentShader=defaultFragmentShaderWithTexture,
                                    camera=view.camera )

mesh._light._pos = (0.5,-0.1,0.5)

# Load mesh to the GPU
view.add(mesh)

# Render
plt.show()
app.run()