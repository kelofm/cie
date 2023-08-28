# --- Python imports ---
import numpy as np
from warnings import warn
from PIL import Image

# --- Vispy imports ---
from vispy.gloo import TextureAtlas

# --- Internal imports ---
from . import textureFolderPath, loadTexture, emptyTexture
from . import incompleteAggregatedVertexShader, incompleteAggregatedFragmentShader

# -----------------------------------------------------------
def checkID(obj,id):
    '''
    Check if the input ID is a valid, which means it is:
        - an integer
        - non-negative
    Issue a warning if the ID is unset (None)
    Raise an exception if the ID does not satisfy ID conditions and is not None
    '''
    if id is None:
        warn(   'Unset ID of object ' + obj.__repr__(),
                category=RuntimeWarning )
    elif not isinstance(id,int) or id<0:
        raise ValueError( 'Invalid ID of object ' + obj.__repr__() )

# -----------------------------------------------------------
class AggregatedTexture(TextureAtlas):
    def __init__(self, **kwargs):
        TextureAtlas.__init__(self, **kwargs)

        # Initialize map:
        #   Stores the base location, width and height of every texture
        #   in the TextureAtlas, accessible by their texture IDs.
        self._map       = {}
        self._IDMap     = np.empty(0)
        self._textures  = []

        # Material constants
        self._ambientMaterialConstants  = []
        self._diffuseMaterialConstants  = []
        self._specularMaterialConstants = []

    def registerTexture(self,   texture,
                                textureID,
                                ambientMaterialConstant=0.5,
                                diffuseMaterialConstant=0.5,
                                specularMaterialConstant=0.5):
        # Check arguments
        if not isinstance(texture, np.ndarray) or not isinstance(texture[0][0][0],np.float32):
            raise TypeError('Textures must be of np.ndarray type with np.float32 type components!')
        
        checkID(texture, textureID)

        # Allocate and check whether it was successful
        bounds = self.get_free_region( texture.shape[0], texture.shape[1] )
        if bounds is None:
            raise RuntimeError('Texture allocation failed!')
        
        # Insert texture
        self._textures.append(texture)
        #self.set_data(np.transpose(texture,axes=(1,0,2)), offset=bounds[0:2], copy=False)

        # Insert material constants
        self._ambientMaterialConstants.append( ambientMaterialConstant )
        self._diffuseMaterialConstants.append( diffuseMaterialConstant )
        self._specularMaterialConstants.append( specularMaterialConstant )

        # Register ID
        if textureID not in list(self._map.keys()):
            self._map[textureID] = bounds
        else:
            raise RuntimeError('The following texture ID is already in use: ' + str(textureID))


    def generateIDMap(self):
        '''
        Create a map that maps the keys of self._map to a 
        contigous sequence 0..n-1
        '''
        if not self._map:
            raise RuntimeError('Cannot generate ID map from empty dict!')

        keys = list( self._map.keys() )
        if max(keys) > len(keys)-1:
            raise RuntimeError('Texture ID is higher than the number of textures!')

        self._IDMap = np.empty(len(keys))
        for index, id in enumerate(keys):
            self._IDMap[id] = index


    def replacePlaceholder(self, string, placeholder, replace, maxIterations=100):
        '''
        Replace every instance of <'placeholder'_placeholder> in 'string' with 'replace'
        '''
        # Check arguments
        if not isinstance(string, str) or not isinstance(placeholder, str) or not isinstance(replace, str):
            raise TypeError('GLSL code, placeholders and their replacements should all be strings!')
        
        # Append placeholder
        placeholder = '<' + placeholder + '_placeholder>'

        # Replace loop
        iterCount   = 0
        while True:
            index   = string.find(placeholder)
            if index is -1:
                break
            else:
                string = string[:index] + replace + string[index+len(placeholder):]

            if iterCount > maxIterations:
                raise RuntimeError('Placeholder replacement exceeded iteration limit!')
            else:
                iterCount += 1

        return string
        

    def generateShaderCode(self):
        # Organize indices and material IDs
        self.generateIDMap()
        # Load shaders
        vertexShader    = incompleteAggregatedVertexShader
        fragmentShader  = incompleteAggregatedFragmentShader
        # Define placeholders
        placeholders    = { 'textureSize'       : len(self._IDMap),
                            'textureIndex'      : [],
                            'textureOffsetX'    : [],
                            'textureOffsetY'    : [],
                            'textureWidth'      : [],
                            'textureHeight'     : [],
                            'ambientConstants'  : [],
                            'diffuseConstants'  : [],
                            'specularConstants' : [] }
        # Fill placeholder values
        for index, item in enumerate( list( self._map.values() ) ):
            placeholders['textureIndex'].append( int( self._IDMap[index] ) )
            placeholders['textureOffsetX'].append( item[0] )
            placeholders['textureOffsetY'].append( item[1] )
            placeholders['textureWidth'].append( item[2] )
            placeholders['textureHeight'].append( item[3] )
        for index, (ambient,diffuse,specular) in enumerate(zip( self._ambientMaterialConstants,
                                                                self._diffuseMaterialConstants,
                                                                self._specularMaterialConstants)):
            placeholders['ambientConstants'].append( ambient )
            placeholders['diffuseConstants'].append( diffuse )
            placeholders['specularConstants'].append( specular )
        # Convert to string and strip brackets
        for key, value in list(placeholders.items()):
            if isinstance(value,list) or isinstance(value,tuple):
                strValue    = ''
                specifier   = 'u,' if key.find('Constant') is -1 else ','
                for item in value:
                    strValue = strValue + str(item) + specifier
                placeholders[key] = strValue
                #value = str(value)
                #placeholders[key] = value[1:-1]
            else:
                placeholders[key] = str(value) + ','
            # Remove last ,
            placeholders[key]   = placeholders[key][:-1]
        # Replace placeholders
        for placeholder, value in list( placeholders.items() ):
            vertexShader    = self.replacePlaceholder(vertexShader, placeholder, value)
            fragmentShader  = self.replacePlaceholder(fragmentShader, placeholder, value)
        # Return shader code
        return vertexShader, fragmentShader


    def build(self,write=False):
        '''
        Copy the textures into an aggregated texture and set the data in the wrapped object
        Returns itself
        '''
        # Prealloc
        aggregate   = np.zeros( self._shape, dtype=np.float32 )

        # Fill data
        for texture, bounds in zip(self._textures, list(self._map.values())):
            aggregate[ bounds[0]:bounds[0]+bounds[2], bounds[1]:bounds[1]+bounds[3], : ] = texture

        # Set data
        self.set_data(aggregate, copy=True)

        # Write data (for debugging)
        if write:
            img         = Image.new('RGB', list(aggregate.shape)[0:2])
            for i in range( list(aggregate.shape)[0] ):
                for j in range( list(aggregate.shape)[1] ):
                    img.putpixel( (i,j), tuple( aggregate[i,j,:]*255 ) )
            img.save( textureFolderPath + '/AggregatedTexture.png' )
        
        del self._textures

        return self
