# --- External Imports ---
import numpy
from plotly import graph_objects

# --- STD Imports ---
import pathlib
import sys

samples: list[list[float]] = []
with open(pathlib.Path(__file__).absolute().parent / sys.argv[1], "r") as file:
    for line in file:
        if line:
            data: list[str] = line.split(",")
            samples.append([float(v) for v in data])

if len(samples[0]) == 2:
    from matplotlib import pyplot
    figure = pyplot.figure()
    axes = figure.subplots()
    axes.plot([x for x, _ in samples], [y for _, y in samples], "+-")
    pyplot.show()
elif len(samples[0]) == 3:
    resolution: int = int(numpy.sqrt(len(samples)))
    if resolution * resolution == len(samples):
        grid = numpy.array(samples).reshape((resolution, resolution, 3))
        figure = graph_objects.Figure(
            data = graph_objects.Surface(
                x = grid[:, :, 0],
                y = grid[:, :, 1],
                z = grid[:, :, 2]))
        figure.show()
    else:
        from matplotlib import pyplot
        figure = pyplot.figure()
        axes.scatter([x for x, _, _ in samples],
                     [y for _, y, _ in samples],
                     [z for _, _, z in samples])

else:
    raise RuntimeError(f"invalid dimension {len(samples[0] - 1)}")

