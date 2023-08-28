# --- Internal Imports ---
import pyfem.pyfem as fem

# --- Python Imports ---
import unittest
import numpy as np


class StaticLinearHeatElement1DLinear(unittest.TestCase):

    def test_integrateStiffness( self ):
        ''' Test stiffness integration for a single element '''
        # Define stiffness matrix
        stiffness = [
            [ 0.0, 0.0, 0.0 ],
            [ 0.0, 0.0, 0.0 ],
            [ 0.0, 0.0, 0.0 ]
        ]

        # Define updater
        def updateStiffness( i, j, value ):
            stiffness[i][j] += value

        # Instantiate an element
        materialParameter   = fem.makeConstantFunction(2.0)
        loadFunction        = fem.makeConstantFunction(1.0)
        domain              = [ 0.0, 1.0 ]
        DoFs                = [ 0, 2 ]
        element             = fem.StaticLinearHeatElement1DLinear(  materialParameter.get(),
                                                                    loadFunction.get(),
                                                                    domain,
                                                                    DoFs )

        # Test integration
        element.integrateStiffness(updateStiffness)
        for i,row in enumerate(stiffness):
            for j,value in enumerate(row):
                if not (i in DoFs) or not (j in DoFs):
                    self.assertEqual( value, 0.0 )
                elif i==j:
                    self.assertEqual( value, materialParameter(0) )
                else:
                    self.assertEqual( value, -materialParameter(0) )

    
    def test_integrateStiffness_list( self ):
        ''' Test parallel stiffness integration for a list of elements '''
        # Define elements
        numberOfElements    = 150
        elements            = fem.StaticLinearHeatElement1DLinear_List()
        materialParameter   = fem.makeConstantFunction(2.0)
        loadFunction        = fem.makeConstantFunction(1.0)
        loads               = [ loadFunction for i in range(numberOfElements) ]
        domains             = [ [float(i),float(i)+1.0] for i in range(numberOfElements) ]
        DoFs                = [ [i,i+1] for i in range(numberOfElements) ]
        for elementIndex in range(numberOfElements):
            elements.emplace_back(  materialParameter.get(),
                                    loads[elementIndex].get(),
                                    domains[elementIndex],
                                    DoFs[elementIndex] )

        # Define stiffness matrix
        stiffness           = np.zeros( (numberOfElements+1,numberOfElements+1) )

        # Define updater
        def updateStiffness( i, j, value ):
            stiffness[i][j] += value
        print("List integrate")
        fem.integrateStiffness( elements, updateStiffness )
        print("List integrate end")
        
        for i,row in enumerate(stiffness):
            for j,value in enumerate(row):
                if i == j:                                  # <-- main diagonal
                    if i!=DoFs[0][0] and i!=DoFs[-1][-1]:       # <-- internal DoF
                        self.assertEqual( value, 2.0*materialParameter )
                    else:                                       # <-- boundary DoF
                        self.assertEqual( value, materialParameter )
                elif abs(i-j) == 1:                         # <-- first band
                    self.assertEqual( value, -materialParameter )
                else:                                       # <-- non-DoF
                    self.assertEqual( value, 0.0 )
        


if __name__ == "__main__":
    unittest.main()