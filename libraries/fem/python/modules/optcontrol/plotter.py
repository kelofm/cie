
from autograd import numpy as np
import matplotlib.pyplot as plt
import time

class plotter( object ):
  
  plotHandles = []

  def __init__(self, labels, axisBounds):
    plt.ion()
    self.fig, self.ax = plt.subplots()
    self.plotHandles = []
    for l in labels:
      p, = self.ax.plot([],[], label=l)
      self.plotHandles.append( p )
    self.ax.set_title("")
    self.ax.axis(axisBounds)
    self.ax.set_xlabel("coordinate/length")
    self.ax.set_ylabel("value")
    self.ax.grid(True)
    plt.legend(handles=self.plotHandles, loc=0)

  def plot(self, n, p, h, basisFunctions, u, locationMaps, plotHandleIndex):

    nResultPointsPerElement = 5*(p+1)
    xPlot = np.zeros(nResultPointsPerElement*n)
    yPlot = np.zeros(nResultPointsPerElement*n)

    for iElement in xrange(0, n):
      x0 = iElement * h
      lm = locationMaps[iElement]
      dofs = u[lm]
      for iResultPoint in xrange( 0, nResultPointsPerElement ):
        r = 2.0 * iResultPoint / ( nResultPointsPerElement - 1.0 ) - 1.0
        x = x0 + ( r + 1.0 ) / 2.0 * h
        xPlot[nResultPointsPerElement*iElement + iResultPoint] = x
        N0 = basisFunctions( r, p, 0 )
        # Evaluate solution
        yPlot[nResultPointsPerElement*iElement + iResultPoint] =  np.dot( N0, dofs ) 
    
    self.plotHandles[plotHandleIndex].set_xdata(xPlot)
    self.plotHandles[plotHandleIndex].set_ydata(yPlot)
  
  def update(self, title, pausetime):
    self.ax.set_title(title)
    self.fig.canvas.draw()
    self.fig.canvas.flush_events()
    time.sleep(pausetime)
    
  def close(self):
    plt.close(self.fig)