# --- External Imports ---
import numpy

# --- STD Imports ---
import abc


class Partitioner(metaclass = abc.ABCMeta):
    @abc.abstractstaticmethod
    def partition(matrix: numpy.array) -> tuple[numpy.array,numpy.array]:
        pass

    @staticmethod
    def _makeRhs(matrix: numpy.array, lhs: numpy.array) -> numpy.array:
        return lhs - matrix


class JacobiPartitioner(Partitioner):
    @staticmethod
    def partition(matrix: numpy.array) -> tuple[numpy.array,numpy.array]:
        lhs = numpy.diag(numpy.diag(matrix))
        return lhs, Partitioner._makeRhs(matrix, lhs)


class GaussSeidelPartitioner(Partitioner):
    @staticmethod
    def partition(matrix: numpy.array) -> tuple[numpy.array,numpy.array]:
        lhs = numpy.diag(numpy.diag(matrix)) + numpy.tril(matrix, -1)
        return lhs, Partitioner._makeRhs(matrix, lhs)
