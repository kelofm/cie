#!/usr/bin/python

# defines useful things for functors
#
# FunctorBinaryOperation(functor1, functor1, operation):
#   Callable object performing given operation on the result from the two functors.
#
# FunctorUnaryOperation(functor, operation):
#   Callable object performing given operation on the result from the functor.
#
# provided binary operations: DifferenceSquared
# provided unary operations: SquareRoot

import math

def DifferenceSquared( x1, x2 ):
  diff = x1 - x2
  return diff * diff

def SquareRoot( x ):
  return math.sqrt( x )

class FunctorBinaryOperation(object):
  def __init__(self, functor1, functor2, operation ):
    self.f1 = functor1
    self.f2 = functor2
    self.op = operation

  def __call__(self, x):
    return self.op( self.f1(x), self.f2(x) )

class FunctorUnaryOperation(object):
  def __init__(self, functor, operation ):
    self.f = functor
    self.op = operation

  def __call__(self, x):
    return self.op( self.f(x) )
