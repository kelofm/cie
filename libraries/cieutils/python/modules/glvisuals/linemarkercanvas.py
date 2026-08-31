# --- Python imports ---
import numpy as np

# --- Vispy imports ---
from vispy import scene

# --- Internal imports ---
from . import ModularCanvas, EventHandlerClass
from . import MarkerVisual, LineVisual

# -----------------------------------------------------------
class LineMarkerCanvas(ModularCanvas):

    def __init__(self):
        ModularCanvas.__init__(self)
        
        self.unfreeze()

        initPos     = np.empty((1,3),dtype=np.float32)
        self.addVisualObject('markers', MarkerVisual(   pos=initPos,
                                                        face_color='w',
                                                        scaling=True)   )
        self.addVisualObject('lines',   LineVisual( pos=initPos ))

        self.eventClass = LineMarkerCanvasCreateEvents(self)
        self.bindEvents(self.eventClass)

        self.freeze()


    def getPos(self,pos):
        transform   = self.scene.node_transform(self.objects['markers'])
        return transform.map(pos)


    def addPoint(self,point):
        self.objects['markers'].addPoint(point)
        self.objects['lines'].addPoint(point)
        self.updateObservers(None)


    def setPoint(self,index,point):
        self.objects['markers'].setPoint(index,point)
        self.objects['lines'].setPoint(index,point)
        self.updateObservers(None)



# -----------------------------------------------------------
class LineMarkerCanvasCreateEvents(EventHandlerClass):
    def __init__(self,*args,**kwargs):
        EventHandlerClass.__init__(self,*args,**kwargs)

    def on_mouse_press(self,event):
        if event.button==1:
            pos = self.canvas.getPos(event.pos)
            self.canvas.addPoint(pos)
        elif event.button==2:
            self.canvas.bindEvents( LineMarkerCanvasEditEvents() )

    def on_mouse_release(self,event):
        self.canvas.objects['markers'].updateMarkers()

    def on_mouse_move(self,event):
        pass


class LineMarkerCanvasEditEvents(EventHandlerClass):
    def __init__(self,*args,**kwargs):
        EventHandlerClass.__init__(self,*args,**kwargs)

    def on_mouse_press(self,event):
        if event.button==1:
            pos = self.canvas.getPos(event.pos)
            self.canvas.objects['markers'].selected_point, self.canvas.objects['markers'].selected_index = self.canvas.objects['markers'].selectPoint(pos)
        elif event.button==2:
            self.canvas.bindEvents( LineMarkerCanvasCreateEvents() )

    def on_mouse_release(self,event):
        # Clear selection
        self.canvas.objects['markers'].selected_point = None
        self.canvas.objects['markers'].selected_index = -1
        self.canvas.objects['markers'].updateMarkers()

    def on_mouse_move(self,event):
        if event.button == 1:   
            self.on_mouse_drag(event)

    def on_mouse_drag(self,event):
        pos = self.canvas.getPos(event.pos)
        if self.canvas.objects['markers'].selected_point is not None:
            # Snap to grid
            self.canvas.objects['markers'].selected_point = pos
            # Update polygon
            self.canvas.setPoint(self.canvas.objects['markers'].selected_index,self.canvas.objects['markers'].selected_point)
            self.canvas.objects['markers'].updateMarkers(self.canvas.objects['markers'].selected_index)