# --- Python imports ---
import numpy as np
from warnings import warn

# --- Vispy imports ---
from vispy.geometry import MeshData

# --- Internal imports ---
from . import checkID

# -----------------------------------------------------
class MeshNode(MeshData):
    def __init__(self, parent=None, children=None, objectID=0, materialID=0, textureCoordinates=None, **kwargs):
        MeshData.__init__(self,**kwargs)

        # Set parent (None=root)
        self._parent    = parent

        # Set children
        if children != None:
            if isinstance(children,list):
                self._children  = children
            else:
                self._children  = [children]
        else:
            self._children  = []

        # Set IDs
        self._objectID      = objectID
        self._materialID    = materialID

        checkID(self,self._objectID)
        checkID(self,self._materialID)

        # Set empty array for unset data, ndarray for set data
        if self._vertices is None:
            self._vertices          = np.empty(0, dtype=np.float32)
            warn('Unset vertices for object ' + self.__repr__(), category=RuntimeWarning)
        elif isinstance( self._vertices, list ) or isinstance( self._vertices, tuple ):
            self._vertices          = np.asarray( self._vertices, dtype=np.float32 )
        elif not isinstance( self._vertices, np.ndarray ):
            raise ValueError('Invalid container type!')

        if self._faces is None:
            self._faces             = np.empty(0, dtype=np.uint32)
            warn('Unset faces for object ' + self.__repr__(), category=RuntimeWarning)
        elif isinstance( self._faces, list ) or isinstance( self._faces, tuple ):
            self._faces             = np.asarray( self._faces, dtype=np.uint32 )
        elif not isinstance( self._faces, np.ndarray ):
            raise ValueError('Invalid container type!')

        if textureCoordinates is None:
            self._textureCoordinates    = np.zeros( (len(self._vertices),2) , dtype=np.float32)
            warn('Unset texture coordinates for object ' + self.__repr__(), category=RuntimeWarning)
        elif isinstance(textureCoordinates, list) or isinstance(textureCoordinates, tuple):
            self._textureCoordinates    = np.asarray( textureCoordinates, dtype=np.float32 )
        elif isinstance(textureCoordinates, np.ndarray):
            self._textureCoordinates    = textureCoordinates
        else:
            raise ValueError('Invalid container type!')


    def parent(self):
        return self._parent


    def setParent(self, parent):
        self._parent = parent


    def child(self,index):
        return self._children[index]


    def children(self):
        return self._children


    def addChild(self,child):
        self._children.append(child)
        child.setParent(self)


    def get_texture_coordinates(self, indexed=None):
        if indexed is None:
            return self._textureCoordinates
        elif indexed is 'faces':
            return self._textureCoordinates[ self._faces ]
        else:
            Exception("Invalid indexing mode. Accepts: None, 'faces'")


    def computeCompiledMeshSize(self):
        '''
        Recursively compute the number of vertices in the entire mesh tree
        '''
        # Prealloc
        sizes   = { 'faces'         : 0,
                    'vertices'      : 0,
                    'nodes'         : 0}

        # Accumulate children sizes
        if self._children:
            for child in self._children:
                childSizes  = child.computeCompiledMeshSize()
                for name, value in childSizes.items():
                    sizes[name] += value

        # Add self sizes
        sizes['faces']          += len(self._faces)
        sizes['vertices']       += len(self._faces)*3
        sizes['nodes']          += 1

        return sizes


    def getCompiledMesh(self, destinationDict=None, faceIndex = 0):
        '''
        Get all vertices, normals, texture coordinates and IDs in the mesh
        '''
        # Create destination
        returnFlag  = False
        if destinationDict is None:
            returnFlag      = True
            sizes           = self.computeCompiledMeshSize()
            destinationDict = { 'vertices'              : np.empty( (sizes['faces'],3,3),  dtype=np.float32),
                                'normals'               : np.empty( (sizes['faces'],3,3),  dtype=np.float32),
                                'textureCoordinates'    : np.empty( (sizes['faces'],3,2),  dtype=np.float32),
                                'objectIDs'             : np.empty( (sizes['faces'],3,1),  dtype=np.float32 ),
                                'materialIDs'           : np.empty( (sizes['faces'],3,1),  dtype=np.float32 )}

        # Write to container - vertices, normals, IDs
        for i, (vertices, normals, textureCoordinates) in enumerate( zip(   self.get_vertices(indexed='faces'),
                                                                            self.get_vertex_normals(indexed='faces'),
                                                                            self._textureCoordinates[self.get_faces()] ) ):
            index   = faceIndex + i
            # Write Data
            destinationDict['vertices'][index]                      = vertices
            destinationDict['normals'][index]                       = normals
            destinationDict['textureCoordinates'][index]            = textureCoordinates
            # Write IDs
            destinationDict['objectIDs'][index]                     = self._objectID * np.ones( (3,1), dtype=np.float32 )
            destinationDict['materialIDs'][index]                   = self._materialID * np.ones( (3,1), dtype=np.float32 )
            '''
            for j, (vertex, normal, coordinates) in enumerate( zip(vertices,normals, textureCoordinates) ):
                index   = vertexIndex + 3*i + j
                # Write data
                destinationDict['vertices'][index]                  = vertex
                destinationDict['normals'][index]                   = normal
                destinationDict['textureCoordinates'][index]        = coordinates
                # Write IDs
                destinationDict['objectIDs'][index]                 = np.uint32(self._objectID)
                destinationDict['materialIDs'][index]               = np.uint32(self._materialID)
        '''
        # Recursive call on the children
        faceIndex   += len(self._faces)
        for child in self._children:
            faceIndex = child.getCompiledMesh(  destinationDict=destinationDict,
                                                faceIndex=faceIndex )

        # Return container if requested
        if returnFlag:
            return destinationDict
        else:
            return faceIndex



# -----------------------------------------------------
class TexturedMeshNode(MeshNode):
    def __init__(self,texture,**kwargs):
        MeshNode.__init__(self,**kwargs)
        # Set texture
        if isinstance( texture, np.ndarray ):
            pass