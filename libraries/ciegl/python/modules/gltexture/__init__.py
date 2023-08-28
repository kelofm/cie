try:
    from cie.__ciegl import *
except:
    pass

from .textureimport import textureFolderPath
from .textureimport import loadTexture, emptyTexture, uniformTexture

from .vertexshaders import defaultVertexShaderWithTexture, incompleteAggregatedVertexShader
from .fragmentshaders import defaultFragmentShaderWithTexture, incompleteAggregatedFragmentShader

from .texture import AggregatedTexture