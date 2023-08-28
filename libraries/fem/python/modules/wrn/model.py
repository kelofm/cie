# --- External Imports ---
import numpy
import scipy

# --- Internal Imports ---
from .interpolation import interpolateFill

# --- STD Imports ---
import typing


class Solution:
    def __init__(self, timeSamples: numpy.array, numberOfVariables: int, stepFrequencies: list[int] = None):
        self.__timeSamples = timeSamples.copy()
        self.__values = numpy.zeros((len(timeSamples), 3, numberOfVariables)) # [i_step, i_derivative, i_variable]
        self.__stepFrequencies = stepFrequencies.copy() if stepFrequencies else [1 for dummy in range(numberOfVariables)]

    @property
    def timeSamples(self) -> numpy.array:
        return self.__timeSamples

    def getTimeSamples(self, i_variable: int) -> numpy.array:
        return self.__timeSamples[::self.getStepFrequency(i_variable)]

    @property
    def values(self) -> numpy.array:
        return self.__values

    @property
    def stepFrequencies(self) -> list[int]:
        return self.__stepFrequencies

    @property
    def numberOfVariables(self) -> int:
        return len(self.__stepFrequencies)

    @property
    def numberOfSteps(self) -> int:
        """Number of steps for the fastest variable."""
        return len(self.timeSamples)

    @property
    def stepSizeBase(self) -> float:
        return self.__timeSamples[1] - self.__timeSamples[0]

    def getStepSize(self, i_variable: int) -> float:
        return self.stepSizeBase * self.getStepFrequency(i_variable)

    def getStepFrequency(self, i_variable: int) -> int:
        return self.__stepFrequencies[i_variable]

    def getNumberOfSteps(self, i_variable: int) -> int:
        return (len(self.timeSamples) - 1) // self.getStepFrequency(i_variable) + 1

    def get(self, i_variable: int, derivative: int = 0, dense: bool = False) -> numpy.array:
        return self.__values[::(1 if dense else self.getStepFrequency(i_variable)), derivative, i_variable]

    def getVariables(self, i_step: int, derivative: int = 0) -> numpy.array:
        return self.values[i_step, derivative, :]

    def set(self, timeHistory: numpy.array, i_variable: int, derivative: int = 0, dense = False) -> None:
        if dense:
            self.__values[:, derivative, i_variable] = timeHistory
        else:
            self.__values[::self.getStepFrequency(i_variable), derivative,  i_variable] = timeHistory
            interpolateFill(self.__timeSamples, self.__values[:, derivative, i_variable], self.getStepFrequency(i_variable))

    def copy(self) -> "Solution":
        copy =  Solution(self.timeSamples, self.numberOfVariables, self.stepFrequencies)
        copy.values[:] = self.__values
        return copy


class Model:
    def __init__(self,
                 massMatrix: numpy.array,
                 dampingMatrix: numpy.array,
                 stiffnessMatrix: numpy.array,
                 loadGenerator: typing.Callable,
                 initialState: tuple[list[float],list[float],list[float]]):
        """ @brief Construct a model from its structural components and initial state.
            @param massMatrix
            @param dampingMatrix
            @param stiffnessMatrix
            @param loadGenerator: functor taking a step (int) and time (float) and returning the current load vector.
            @param initialState: [values, derivatives, secondDerivatives].
        """
        self.__massMatrix = massMatrix
        self.__dampingMatrix = dampingMatrix
        self.__stiffnessMatrix = stiffnessMatrix

        self.__loadGenerator: typing.Callable = None
        self.setLoadGenerator(loadGenerator)

        self.__massInverse: numpy.array = None
        self.__dampingInverse: numpy.array = None
        self.__stiffnessInverse: numpy.array = None

        self.__initialState: numpy.array = None
        self.setInitialState(initialState)


    @property
    def massMatrix(self) -> numpy.array:
        return self.__massMatrix

    @property
    def dampingMatrix(self) -> numpy.array:
        return self.__dampingMatrix

    @property
    def stiffnessMatrix(self) -> numpy.array:
        return self.__stiffnessMatrix

    @property
    def massInverse(self) -> numpy.array:
        if self.__massInverse is None:
            self.__massInverse = scipy.linalg.inv(self.massMatrix)
        return self.__massInverse

    @property
    def dampingInverse(self) -> numpy.array:
        if self.__dampingInverse is None:
            self.__dampingInverse = scipy.linalg.inv(self.__dampingMatrix)
        return self.__dampingInverse

    @property
    def stiffnessInverse(self) -> numpy.array:
        if self.__stiffnessInverse is None:
            self.__stiffnessInverse = scipy.linalg.inv(self.stiffnessMatrix)
        return self.__stiffnessInverse

    @property
    def loadGenerator(self) -> numpy.array:
        return self.__loadGenerator

    @loadGenerator.setter
    def loadGenerator(self, generator: typing.Callable) -> None:
        self.__loadGenerator = generator

    @property
    def initialState(self) -> numpy.array:
        return self.__initialState

    @property
    def numberOfVariables(self) -> int:
        return self.massMatrix.shape[0]

    def setLoadGenerator(self, loadGenerator: typing.Callable) -> None:
        self.__loadGenerator = loadGenerator

    def setInitialState(self, state: tuple[list[float],list[float],list[float]]) -> None:
        # Append initial conditions with initial acceleration
        initialAccelerations = self.massInverse.dot(self.makeLoadVector(0, 0.0) - self.stiffnessMatrix.dot(state[0]) - self.dampingMatrix.dot(state[1]))
        state = list(state) + [initialAccelerations]
        self.__initialState = numpy.array(state)

    def makeLoadVector(self, step: int, time: float) -> numpy.array:
        return self.__loadGenerator(step, time)

    def applyInitialConditions(self, values: numpy.array, derivatives: numpy.array, secondDerivatives: numpy.array) -> None:
        values[0,:] = self.__initialState[0,:]
        derivatives[0,:] = self.__initialState[1,:]
        secondDerivatives[0,:] = self.__initialState[2,:]

    def applyInitialConditionsToSolution(self, solution: Solution) -> None:
        self.applyInitialConditions(solution.values[:,0,:], solution.values[:,1,:], solution.values[:,2,:])
