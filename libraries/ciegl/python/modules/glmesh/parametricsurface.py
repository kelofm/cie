# --- Python imports ---
import numpy as np

# --- Internal imports ---
from . import convertToSurfaceMesh

# -----------------------------------------------------------
class ParametricSurface:
    def __init__(self, vertexFunctor, domain=((0.0,1.0),(0.0,1.0))):        
        self._functor   = vertexFunctor
        self._domain    = domain


    def evaluate(self, u, v, format='grid'):
        # Check arguments
        if min(u)<min(self._domain[0]) or max(u)>max(self._domain[0]):
            raise ValueError('Parameter u out of range')
        if min(v)<min(self._domain[1]) or max(v)>max(self._domain[1]):
            raise ValueError('Parameter v out of range')

        if format is 'grid':
            return np.array(
                    [
                        [ [self._functor[componentIndex](U,V) for U in u] for V in v ] for componentIndex in (0,1,2)
                    ]
                    ,dtype=np.float32)
        elif format is 'list':
            return convertToSurfaceMesh(
                        np.array(
                            [
                                [ [self._functor[componentIndex](U,V) for U in u] for V in v ] for componentIndex in (0,1,2)
                            ]
                            ,dtype=np.float32)
                    )
        else:
            raise ValueError('Invalid output format')


    def u(self,nSamples=10):
        return np.linspace( self._domain[0][0], self._domain[0][1], num=nSamples )


    def v(self,nSamples=10):
        return np.linspace( self._domain[1][0], self._domain[1][1], num=nSamples )


# -----------------------------------------------------------
Torus       = lambda radius1, radius2, offset : ParametricSurface(    (   
                                                    lambda u,v : np.cos(u)*(radius1+radius2*np.cos(v)) + offset[0],
                                                    lambda u,v : np.sin(u)*(radius1+radius2*np.cos(v)) + offset[1],
                                                    lambda u,v : radius2*np.sin(v) + offset[2] ),
                                                    domain=((0.0,2*np.pi),(0.0,-2*np.pi))                           )


Sphere      = lambda radius, offset : ParametricSurface(    (  lambda u,v : radius * np.cos(u) * np.cos(v) + offset[0],
                                                               lambda u,v : radius * np.sin(u) * np.cos(v) + offset[1],
                                                               lambda u,v : radius * np.sin(v) + offset[2]  ),
                                                            domain=((0.0, 2*np.pi), (np.pi/2.0,-np.pi/2.0))              )