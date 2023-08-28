# --- External Imports ---
import numpy

# --- Internal Imports ---
from .partitioning import Partitioner, JacobiPartitioner, GaussSeidelPartitioner
from .time_integration import TimeIntegrator, NewmarkTimeIntegrator
from .model import Model, Solution

# --- STD Imports ---
import abc
import functools


class WaveformRelaxation(metaclass = abc.ABCMeta):
    def __init__(self, model: Model):
        self.__model = model
        self.__partitionedModel: Model = None

    @abc.abstractproperty
    def partitioner(self) -> Partitioner:
        pass

    @abc.abstractproperty
    def solution(self) -> Solution:
        pass

    @abc.abstractmethod
    def getTimeIntegrator(self, *args, **kwargs) -> TimeIntegrator:
        pass

    @abc.abstractmethod
    def relax(self) -> float:
        """@brief Perform a relax step and return the solution update's norm."""
        pass

    def solve(self, updateTolerance: float = 1e-10, maxSteps: int = int(1e6), hook = lambda this, **kwargs: False) -> bool:
        norm = float("inf")
        i_relax = 0
        while updateTolerance < norm and i_relax < maxSteps:
            norm = self.relax()
            hook(norm = norm, i_relax = i_relax)
            i_relax += 1

    @property
    def model(self) -> Model:
        return self.__model

    @property
    def partitionedModel(self) -> Model:
        if not self.__partitionedModel:
            massLhs, massRhs = self.partitioner.partition(self.model.massMatrix)
            dampingLhs, dampingRhs = self.partitioner.partition(self.model.dampingMatrix)
            stiffnessLhs, stiffnessRhs = self.partitioner.partition(self.model.stiffnessMatrix)
            loadGenerator = lambda step, time: self.model.loadGenerator(step, time) \
                                               + massRhs.dot(self.solution.getVariables(step, derivative = 2)) \
                                               + dampingRhs.dot(self.solution.getVariables(step, derivative = 1)) \
                                               + stiffnessRhs.dot(self.solution.getVariables(step))
            self.__partitionedModel = Model(massLhs, dampingLhs, stiffnessLhs, loadGenerator, self.model.initialState)
        return self.__partitionedModel


class HomogeneousRelaxation(WaveformRelaxation):
    def __init__(self,
                 model: Model,
                 solution: Solution,
                 partitioner: Partitioner = GaussSeidelPartitioner()):
        super().__init__(model)
        self.__partitioner = partitioner
        self.__solution = solution

        i_minStepFrequency = numpy.argmin(self.solution.stepFrequencies)
        self.__timeIntegrator = NewmarkTimeIntegrator(self.solution.getStepSize(i_minStepFrequency),
                                                      self.solution.getNumberOfSteps(i_minStepFrequency),
                                                      begin = self.solution.timeSamples[0])

    @property
    def solution(self) -> Solution:
        return self.__solution

    @property
    def partitioner(self) -> Partitioner:
        return self.__partitioner

    def getTimeIntegrator(self) -> TimeIntegrator:
        return self.__timeIntegrator

    def relax(self) -> float:
        displacements, velocities, accelerations = self.getTimeIntegrator().integrate(self.partitionedModel)
        norm = numpy.max(numpy.abs(displacements - self.solution.values[:,0,:]))
        self.solution.values[:,0,:] = displacements
        self.solution.values[:,1,:] = velocities
        self.solution.values[:,2,:] = accelerations
        return norm


class HeterogeneousJacobiRelaxation(WaveformRelaxation):
    def __init__(self,
                 model: Model,
                 solution: Solution,
                 timeIntegrators = None):
        super().__init__(model)
        self.__solution = solution

        if timeIntegrators:
            self.__timeIntegrators = timeIntegrators
        else:
            self.__timeIntegrators = []
            for i_variable in range(self.solution.numberOfVariables):
                stepSize = self.solution.getStepSize(i_variable)
                numberOfSteps = self.solution.getNumberOfSteps(i_variable)
                self.__timeIntegrators.append(NewmarkTimeIntegrator(stepSize, numberOfSteps, begin = self.__solution.timeSamples[0]))

    @property
    def solution(self) -> Solution:
        return self.__solution

    @property
    def partitioner(self):
        return JacobiPartitioner()

    def getTimeIntegrator(self, i_variable: int) -> TimeIntegrator:
        return self.__timeIntegrators[i_variable]

    def relax(self) -> float:
        norm = 0

        @functools.cache
        def forces(i_variable: int, step: int, time: float) -> numpy.array:
            return self.partitionedModel.makeLoadVector(self.solution.getStepFrequency(i_variable) * step, time)

        for i_variable in range(self.model.numberOfVariables):
            localMass = numpy.array([[self.partitionedModel.massMatrix[i_variable,i_variable]]])
            localDamping = numpy.array([[self.partitionedModel.dampingMatrix[i_variable,i_variable]]])
            localStiffness = numpy.array([[self.partitionedModel.stiffnessMatrix[i_variable,i_variable]]])
            localForces = lambda step, time: numpy.array([forces(i_variable, step, time)[i_variable]])
            localInitialState = [[derivative[i_variable]] for derivative in self.partitionedModel.initialState]
            localModel = Model(localMass, localDamping, localStiffness, localForces, localInitialState)

            localSolution, localVelocity, localAcceleration = self.getTimeIntegrator(i_variable).integrate(localModel)
            norm = max([norm, numpy.max(numpy.abs(numpy.ravel(localSolution) - self.solution.get(i_variable)))])
            self.solution.set(numpy.ravel(localSolution), i_variable)
            self.solution.set(numpy.ravel(localVelocity), i_variable, derivative = 1)
            self.solution.set(numpy.ravel(localAcceleration), i_variable, derivative = 2)
        return norm
