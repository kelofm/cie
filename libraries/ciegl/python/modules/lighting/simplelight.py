import numpy as np

class SimpleLight:
    # Simple spot light
    def __init__(self,  parent=None, 
                            pos=[0.0,0.0,0.0], 
                            color=[1.0,1.0,1.0],
                            ambient=[0.2,0.2,0.2] ):

        self._parent        = parent
        self._pos           = np.array( pos, dtype=np.float32 )
        self._color         = np.array( color, dtype=np.float32 )
        self._ambient       = np.array( ambient, dtype=np.float32 )


    def update(self, 
                    pos=None, 
                    color=None,
                    ambient=None ):
        if pos is not None:
            self._pos       = np.array( pos, dtype=np.float32 )
        if color is not None:
            self._color     = np.array( color, dtype=np.float32 )
        if ambient is not None:
            self._ambient   = np.array( ambient, dtype=np.float32 )
        if self._parent is not None:
            self._parent.updateLight()



    def setParent(self,parent):
        self._parent = parent