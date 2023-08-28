# --- Python Imports ---
from matplotlib import pyplot as plt
from matplotlib import animation as anim
import numpy as np

# --- Internal Imports ---


# ---------------------------------------------------------
def animateTimeSeries(  time,
                        positions,
                        timeSeriesSolution,
                        model,
                        speed=0.1,
                        repeat=True,
                        ylim=(-1.0,1.0),
                        figure=None,
                        axis=None ):
    # Pack into list in only one time series is supplied
    if len( timeSeriesSolution.shape ) < 3:
        timeSeriesSolution = np.asarray([timeSeriesSolution])

    values  = [model.sample( timeSeries[0], positions ) for timeSeries in timeSeriesSolution]
    blit    = True

    # Create figure and axis if not given
    if figure is None or axis is None:
        figure, axis    = plt.subplots()

    # Initialize lines
    lines       = [ None for i in range(len(timeSeriesSolution)) ]
    for index, valueSet in enumerate( values ):
        lines[index], = axis.plot( positions, valueSet )
    timeLabel   = axis.text(    0.01,
                                0.95,
                                "t = %.3f" % time[0],
                                fontdict={  "family"    : "serif",
                                            "alpha"     : 0.2,
                                            "weight"    : "bold"} )
    axis.set_ylim(ylim)

    # Set callbacks
    def frameGenerator():
        for index, t in enumerate(time):
            yield ( t, [model.sample( timeSeries[index], positions ) for timeSeries in timeSeriesSolution] )

    def initAnim():
        for line in lines:
            line.set_ydata( [np.nan] * len(positions) )
        timeLabel.set_text( "t = %.3f" % time[0] )
        return (*lines, timeLabel)

    def stepAnim( args ):
        # Update time
        timeLabel.set_text( "t = %.3f" % args[0] )

        # Update plot
        for index, line in enumerate( lines ):
            line.set_ydata( args[1][index] )
        return (*lines, timeLabel)

    # Create animation
    animation = anim.FuncAnimation( figure,
                                    stepAnim,
                                    init_func=initAnim,
                                    frames=frameGenerator,
                                    interval=1.0/speed,
                                    blit=blit,
                                    save_count=50,
                                    repeat_delay=1000,
                                    repeat=repeat)
    return animation




class ConvergencePlot:
    '''
    Convergence plot that can be updated dynamically
    '''
    def __init__( self, maxPoints=50, figure=None, axes=None, yLabel="value [-]" ):
        if figure is None and axes is None:
            figure  = plt.figure()
            axes    = figure.add_subplot( 1,1,1 )
        elif (figure is None) is not (axes is None):
            raise RuntimeError( "Either BOTH figure AND axes must be given, or NEITHER" )

        self.maxPoints  = maxPoints
        self.reset()

        self.figure     = figure
        self.axes       = axes
        self.line,      = plt.plot( self._x, self._y, "--." )

        self.figure.suptitle( "Convergence Plot" )
        self.axes.set_yscale( "log" )
        self.axes.set_xlabel( "iteration#" )
        self.axes.set_ylabel( yLabel )

        self.start()


    def reset( self ):
        self._x         = []
        self._y         = []
        self._yMax      = 0
        self._yMin      = None
        self._counter   = 0


    def add( self, value ):
        self._counter += 1
        if self._yMax < value:
            self._yMax = value

        if value !=0.0 and (self._yMin is None or self._yMin > value):
            self._yMin = value

        if self._counter <= self.maxPoints:
            self._x.append( self._counter )
            self._y.append( value )
        else:
            self._x[0] += 1
            for k in range(1, self.maxPoints):
                self._x[k]      += 1
                self._y[k-1]    = self._y[k]
            self._y[-1] = value



    def __call__( self, value ):
        self.add( value )
        self.line.set_xdata( self._x )
        self.line.set_ydata( self._y )
        self.axes.set_xlim( self._x[0], self._x[-1] )
        self.axes.set_ylim( self._yMin, self._yMax )
        self.figure.canvas.draw()
        plt.pause(1e-3)


    def start( self ):
        plt.ion()

    def finish( self ):
        plt.ioff()


    def close( self ):
        self.finish()
        plt.close( self.figure )