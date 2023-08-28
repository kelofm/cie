# --- External Imports ---
import numpy

# --- STL Imports ---
import json


if __name__ == "__main__":

    lhs = lambda x, y, z: -3.0 - x + y -5.0*z
    rhs = lambda x, y, z: numpy.array([ lhs(x, y, z), -3.0 -5.0*x + y - z ])

    dLhs = lambda x, y, z: numpy.array([ -1.0, 1.0, -5.0 ])
    dRhs = lambda x, y, z: numpy.array([ [-1.0, 1.0, -5.0], [-5.0, 1.0, -1.0] ])

    product = lambda x, y, z: lhs(x, y, z) * rhs(x, y, z)
    dProduct = lambda x, y, z: lhs(x, y, z) * dRhs(x, y, z) + numpy.outer( rhs(x, y, z), dLhs(x, y, z) )

    results = []

    for i in range( -2, 3 ):
        for j in range( -2, 3 ):
            for k in range( -2, 3 ):
                results.append({ 
                    "point" : [i, j, k],
                    "value" : list(product(i, j, k)),
                    "derivative" : list(numpy.ravel(dProduct(i, j, k)))
                })

    with open( "scalarVectorProductReference.json", "w" ) as file:
        json.dump( {"results" : results}, file, indent="    " )