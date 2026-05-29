# -----------------------------------------------------------
class Observer:
    def __init__(self,subject):
        self.subject = subject

    def attach(self,subject):
        self.subject = subject

    def detach(self):
        self.subject = None

    def onSubjectChange(self,event):
        pass