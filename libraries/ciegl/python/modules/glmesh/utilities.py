# --- Python imports ---
from warnings import warn

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