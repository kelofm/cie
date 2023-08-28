# --- External Imports ---
import numpy

# --- Internal Imports ---
from .model import Solution, Model


class WindowSequence:
    def __init__(self, model: Model, completeSolution: Solution, numberOfWindows: int):
        self.__model = model
        self.__completeSolution = completeSolution
        self.__numberOfWindows = numberOfWindows

        self.__currentSlice: numpy.s_ = None
        self.__currentSolution: Solution = None
        self.__windowIndexGenerator = iter(range(numberOfWindows))

        # Apply time offset
        self.__stepOffset = 0
        self.__defaultLoadGenerator = self.__model.loadGenerator
        self.__model.setLoadGenerator(self.__offsetLoadGenerator)

    @property
    def solution(self) -> Solution:
        return self.__completeSolution

    @property
    def windowSolution(self) -> Solution:
        return self.__currentSolution

    @property
    def model(self) -> Model:
        return self.__model

    def updateSolution(self) -> None:
        self.solution.values[self.__currentSlice,:,:] = self.__currentSolution.values

    def __offsetLoadGenerator(self, step: int, time: float) -> numpy.array:
        return self.__defaultLoadGenerator(step - self.__stepOffset, time)

    def __getWindowBegin(self, i_window: int) -> int:
        i_maxStepFrequency = numpy.argmax(self.solution.stepFrequencies)
        slowIndex = (i_window * self.solution.getNumberOfSteps(i_maxStepFrequency) // self.__numberOfWindows - 1) if i_window else 0
        return slowIndex * self.solution.getStepFrequency(i_maxStepFrequency)

    def __getWindowEnd(self, i_window: int) -> int:
        """Last window may be longer than the others."""
        if i_window == self.__numberOfWindows - 1:
            return self.solution.numberOfSteps
        else:
            return self.__getWindowBegin(i_window + 1) + 1
            #i_maxStepFrequency = numpy.argmax(self.solution.stepFrequencies)
            #slowIndex = ((i_window + 1) * self.solution.getNumberOfSteps // self.__numberOfWindows)
            #return slowIndex * self.solution.getStepFrequency(i_maxStepFrequency)

    def __nextWindow(self) -> None:
        isFirstWindow = self.__currentSolution == None

        # Negate current offset
        if not isFirstWindow:
            self.updateSolution()
            stepOffset = self.__currentSolution.numberOfSteps - 1

            # Transfer the end state from the current solution
            # to the initial state of the next one.
            i_lastStep = self.__currentSolution.numberOfSteps - 1
            nextInitialState = list(self.__currentSolution.getVariables(i_lastStep, derivative) for derivative in range(3))
            self.__model.setInitialState(nextInitialState)

        # Slice the complete solution for the next window
        i_window = next(self.__windowIndexGenerator)
        self.__currentSlice = numpy.s_[self.__getWindowBegin(i_window):self.__getWindowEnd(i_window)]
        self.__currentSolution = Solution(self.solution.timeSamples[self.__currentSlice],
                                          self.solution.numberOfVariables,
                                          self.solution.stepFrequencies)

        # Apply next offset
        if not isFirstWindow:
            self.__stepOffset += stepOffset

    def __iter__(self) -> "WindowSequence":
        return self

    class Window:
        def __init__(self, solution: Solution, model: Model, stepOffset: int):
            self.__solution = solution
            self.__model = model
            self.__stepOffset = stepOffset

        @property
        def solution(self) -> Solution:
            return self.__solution

        @property
        def model(self) -> Model:
            return self.__model

        @property
        def stepOffset(self) -> int:
            return self.__stepOffset

    def __next__(self) -> Window:
        self.__nextWindow()
        return WindowSequence.Window(self.windowSolution, self.model, self.__stepOffset)
