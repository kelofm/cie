# --- External Imports ---
import numpy
import scipy.linalg

# --- Internal Imports ---
from .model import Model

# --- STD Imports ---
import abc


class TimeIntegrator(metaclass = abc.ABCMeta):
    @abc.abstractmethod
    def integrate(self, model: Model) -> tuple[numpy.array,numpy.array]:
        pass


class NewmarkTimeIntegrator(TimeIntegrator):
    def __init__(self, stepSize: float, numberOfSteps: int, begin: float = 0.0, beta: float = 0.25, gamma: float = 0.5):
        super().__init__()
        self.__begin = begin
        self.__stepSize = stepSize
        self.__stepSize2 = stepSize * stepSize
        self.__numberOfSteps = numberOfSteps
        self.__beta = beta
        self.__gamma = gamma
        if self.__beta != 0.25 or self.__gamma != 0.5:
            raise RuntimeError("Newmark integration is only implemented for beta=0.25 and gamma=0.5")

    def integrate(self, model: Model) -> tuple[numpy.array,numpy.array,numpy.array]:
        # Initialize
        displacements = model.initialState[0,:].copy()
        velocities = model.initialState[1,:].copy()
        accelerations = model.initialState[2,:].copy()

        lhs = model.massMatrix + self.__gamma * self.__stepSize * model.dampingMatrix + self.__beta * self.__stepSize2 * model.stiffnessMatrix

        # Allocate and initialize output
        displacementHistory = numpy.zeros((self.__numberOfSteps, model.numberOfVariables))
        velocityHistory = numpy.zeros(displacementHistory.shape)
        accelerationHistory = numpy.zeros(displacementHistory.shape)
        model.applyInitialConditions(displacementHistory, velocityHistory, accelerationHistory)

        for i_step in range(1, self.__numberOfSteps):
            time = self.__begin + i_step * self.__stepSize

            # Predictor phase
            displacements += self.__stepSize * velocities + self.__stepSize2 / 2.0 * (1.0 - 2.0 * self.__beta) * accelerations
            velocities += self.__stepSize * (1.0 - self.__gamma) * accelerations

            # Solution
            accelerations = scipy.linalg.solve(lhs,
                                               model.makeLoadVector(i_step, time) - model.dampingMatrix.dot(velocities) - model.stiffnessMatrix.dot(displacements))

            # Corrector phase
            displacements += self.__beta * self.__stepSize2 * accelerations
            velocities += self.__gamma * self.__stepSize * accelerations

            # Update output
            displacementHistory[i_step,:] = displacements.copy()
            velocityHistory[i_step,:] = velocities.copy()
            accelerationHistory[i_step,:] = accelerations.copy()

        return displacementHistory, velocityHistory, accelerationHistory
