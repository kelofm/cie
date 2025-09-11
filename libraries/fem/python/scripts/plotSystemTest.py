# --- External Imports ---
import numpy
from plotly import graph_objects

# --- STD Imports ---
import sys
from typing import Optional


if sys.argv[1].endswith("csv"):
    # Settings.
    useMatplotlib: bool = False
    elementPlot: bool   = False

    # Read csv header.
    cellsPerDirection: int
    postprocessResolution: int
    iXColumn: int
    iYColumn: Optional[int]
    iStateColumn: int

    headerLines: list[list[str]]
    with open(sys.argv[1], "r") as file:
        headerLines = [file.readline().strip().split(",") for _ in range(3)]

    cellsPerDirection = int(headerLines[1][headerLines[0].index("cells-per-direction")])
    postprocessResolution = int(headerLines[1][headerLines[0].index("postprocess-resolution")])
    iXColumn = headerLines[2].index("x")
    iYColumn = headerLines[2].index("y")
    iYColumn = iYColumn if iYColumn != -1 else None
    iStateColumn = headerLines[2].index("state")

    dimension: int = 1 if iYColumn is None else 2

    # Load data from the CSV.
    samples: numpy.ndarray

    if dimension == 1:
        samples = numpy.loadtxt(sys.argv[1], delimiter = ",", skiprows = 3)
    elif dimension == 2:
        samples = numpy.loadtxt(
            sys.argv[1],
            delimiter = ",",
            skiprows = 3
        ).reshape((
            numpy.pow(cellsPerDirection, dimension),
            numpy.pow(postprocessResolution, dimension),
            dimension + 1,
        ), order = "C")
    else:
        raise RuntimeError(f"invalid dimension {dimension}")

    # Scatter plot.
    if dimension == 1:
        from matplotlib import pyplot
        figure = pyplot.figure()
        axes = figure.subplots()
        axes.plot([sample[iXColumn] for sample in samples], [sample[iStateColumn] for sample in samples], "+-")
        pyplot.show()
    elif dimension == 2:
        if useMatplotlib:
            from matplotlib import pyplot
            figure = pyplot.figure()
            axes = figure.add_subplot(projection = '3d')
            axes.scatter(samples[:,iXColumn],
                        samples[:,iYColumn],
                        samples[:,iStateColumn])
            pyplot.show()
        else:
            data = []

            if elementPlot and postprocessResolution != 1:
                samples = samples.reshape((numpy.pow(cellsPerDirection, dimension),
                                        postprocessResolution,
                                        postprocessResolution,
                                        dimension + 1))
                data += [
                    graph_objects.Surface(
                        x = samples[iCell, :, :, iXColumn],
                        y = samples[iCell, :, :, iYColumn],
                        z = samples[iCell, :, :, iStateColumn])
                    for iCell in range(samples.shape[0])]
            else:
                samples = samples.reshape((cellsPerDirection * postprocessResolution,
                                    cellsPerDirection * postprocessResolution,
                                    dimension + 1))
                data.append(
                    graph_objects.Scatter3d(
                        x = numpy.ravel(samples[:, :, iXColumn]),
                        y = numpy.ravel(samples[:, :, iYColumn]),
                        z = numpy.ravel(samples[:, :, iStateColumn]),
                        mode = "markers",
                        marker = {
                            "size" : 1,
                            "color" : numpy.ravel(samples[:, :, iStateColumn])}))
            figure = graph_objects.Figure(data = data)
            figure.show()
elif sys.argv[1].endswith("stl"):
    # STL plot.
    import stl
    mesh = stl.mesh.Mesh.from_file(sys.argv[1])

    points = mesh.vectors.reshape(-1,3)  # reshape Nx3x3 to N*3x3 to convert polygons to points
    faces_index = numpy.arange(len(points))  # get indexes of points

    figure = graph_objects.Figure(data = [
        graph_objects.Mesh3d(
            x = points[:,0],  # pass first column of points array
            y = points[:,1],  # .. second column
            z = points[:,2],  # .. third column
            # i, j and k give the vertices of triangles
            i = faces_index[0::3],  # indexes of the points. k::3 means pass every third element
            j = faces_index[1::3],  # starting from k element
            k = faces_index[2::3],  # for example 2::3 in arange would be 2,5,8,...
            opacity = .9
        )
    ])

    figure.show()
