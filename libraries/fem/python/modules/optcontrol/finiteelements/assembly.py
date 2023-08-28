#!/usr/bin/python

from autograd import numpy

def scatterMatrixInto( globalMatrix, elementMatrix, locationMap ):

  size1 = elementMatrix.shape[0]  
  size2 = elementMatrix.shape[0]  

  if size1 != elementMatrix.shape[1]:
    raise ValueError('Element matrix must be square!')

  if size1 != locationMap.shape[0]:
    raise ValueError('Element matrix and location map size do not correspond! Make sure location map has size #dof.')

  globalMatrix[numpy.ix_(locationMap, locationMap)] = globalMatrix[numpy.ix_(locationMap, locationMap)] + elementMatrix;

def scatterVectorInto( globalVector, elementVector, locationMap ):
  size1 = elementVector.shape[0]  

  if size1 != locationMap.shape[0]:
    raise ValueError('Element vector and location map size do not correspond! Make sure location map has size #dof.')

  globalVector[locationMap] = globalVector[locationMap] + elementVector
