# --- Python imports ---
import numpy as np

# --- Vispy imports ---
from vispy import scene

# -----------------------------------------------------------
class ModularCanvas(scene.SceneCanvas):

    def __init__(self):
        scene.SceneCanvas.__init__(self,    keys='interactive',
                                            size=(800, 600))

        self.unfreeze()

        # Initialize backend
        self.view           = self.central_widget.add_view()
        self.view.camera    = scene.PanZoomCamera(  rect=(-100, -100, 200, 200),
                                                    aspect=1.0)

        # Initialize state storage
        self.eventClass     = EventHandlerClass(canvas=self)
        self.objects        = {}
        self.observers      = []

        self.freeze()


    def addVisualObject(self,visualObjectName,visualObject):
        self.objects[visualObjectName]  = visualObject
        self.view.add( self.objects[visualObjectName] )


    def bindEvents(self,eventClass):
        self.eventClass     = eventClass
        eventClass.canvas   = self


    def attachObserver(self,observer):
        self.observers.append(observer)
        observer.attach(self)


    def detachObserver(self,observer):
        for index, peep in enumerate(self.observers):
            if peep is observer:
                peep.detach()
                del self.observers[index]
                break


    def updateObservers(self, event):
        for observer in self.observers:
            observer.onSubjectChange(event)


    # EVENT HANDLER ROUTERS --------------------------------
    def on_mouse_press(self,event):
        self.eventClass.on_mouse_press(event)


    def on_mouse_release(self,event):
        self.eventClass.on_mouse_release(event)


    def on_mouse_move(self,event):
        self.eventClass.on_mouse_move(event)




# EVENT HANDLER CLASSES ---------------------------------
class EventHandlerClass:
    def __init__(self,canvas=None):
        self.canvas = canvas

    def on_mouse_press(self,event):
        pass

    def on_mouse_release(self,event):
        pass

    def on_mouse_move(self,event):
        pass