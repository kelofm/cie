# --- Python Imports ---
import collections
import numpy as np

# ---------------------------------------------------------

class Cache:
    '''
    Wrapper class for a dictionary that is used
    for caching lists of values (or objects)
    '''
    def __init__( self ):
        self.map = dict()


    def get( self, key, hashed=False ):
        if not hashed:
            id = self.hash(key)
        else:
            id = key

        if id in self.map:
            return self.map[id]
        else:
            raise ValueError( "Undefined key in cache: " + str(key) )


    def insert( self, key, value, hashed=False ):
        if not hashed:
            id = self.hash(key)
        else:
            id = key

        if id in self.map:
            return False
        else:
            self.map[id] = value
            return True


    def overwrite( self, key, value, hashed = False ):
        if not hashed:
            id = self.hash(key)
        else:
            id = key

        self.map[id] = value


    def check( self, key, hashed=False ):
        if not hashed:
            id = self.hash(key)
        else:
            id = key

        return id in self.map


    def clear( self ):
        self.map = {}


    def hash( self, item ):
        if isinstance( item, collections.abc.Hashable ):
            return hash(item)
        elif isinstance( item, np.ndarray ):
            return hash( item.data.tobytes() )
        else:
            raise TypeError( "Unhashable item of type: " + str(type(item)) )