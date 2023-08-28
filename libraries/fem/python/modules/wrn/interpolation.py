# --- External Imports ---
import numpy


def interpolateFill(times: numpy.array, values: numpy.array, stepFrequency: int) -> None:
    values[:] = numpy.interp(times, times[::stepFrequency], values[::stepFrequency])
