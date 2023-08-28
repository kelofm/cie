# --- External Imports ---
import numpy
from matplotlib import pyplot

# --- Internal Imports ---
from .Result import Result

# --- STD Imports ---
import typing
import math


colors = {
    "TUMBlue" : [v / 255.0 for v in (0, 82, 147, 0xff)],
    "TUMOrange" : [v / 255.0 for v in (227, 114, 34, 0xff)],
    "warning" : [v / 255.0 for v in (0xff, 0xcc, 0xcc, 0xff)],
    "black" : [0.0, 0.0, 0.0, 1.0],
    "gray" : [0.5, 0.5, 0.5, 1.0]
}


def findMinPrefix(contents: list[Result], quantityName: str):
    minPrefix = (None, None)
    prefixes = set()

    if isinstance(contents, list):
        for results in contents:
            for case in results.cases:
                prefixes.add(case[quantityName].getPrefix())
    elif isinstance(contents, Result):
        for case in contents.cases:
            prefixes.add(case[quantityName].getPrefix())
    else:
        return ValueError(f"Unexpected input type: {type(contents)}")

    for prefix in prefixes:
        exponent = Result.Case.Quantity.getPrefixMap()[prefix]
        if minPrefix[1] == None or minPrefix[2] < exponent:
            minPrefix = (prefix, exponent)

    return minPrefix[0]


def homogenizeUnits(results: Result, quantityName: str, commonPrefix: str = "") -> None:
    if not commonPrefix:
        commonPrefix = findMinPrefix(results, quantityName)

    for case in results.cases:
        case[quantityName].setPrefix(commonPrefix)


def plotQuantity(results: Result,
                 quantityName: str,
                 axes: pyplot.Axes,
                 color: list[float] = colors["TUMBlue"],
                 textColor: list[float] = colors["black"],
                 prefix: str = None,
                 label: str = "",
                 positions: numpy.array = None,
                 barCount: int = 1,
                 barIndex: int = 0,
                 totalWidth: float = 0.8,
                 format = R"{:.2E}") -> typing.Any:
    if not prefix:
        prefix = findMinPrefix(results, quantityName)
    homogenizeUnits(results, quantityName, commonPrefix = prefix)

    if not label:
        label = quantityName

    # Collect values
    values = [case[quantityName].getValue() for case in results.cases]

    if positions is None:
        positions = numpy.arange(len(values))

    # Plot values
    barWidth = totalWidth / barCount
    offsets = positions + (barWidth * (barIndex + 0.5) - totalWidth / 2.0)
    bars = axes.bar(offsets,
                    values,
                    barWidth,
                    label = label,
                    color = color)

    # Draw bar heights
    for rectangle in bars.patches:
        height = rectangle.get_height()
        width = rectangle.get_width()
        position = rectangle.get_x()
        axes.text(
            position + width / 2.0,
            height * ((barIndex + 1) / barCount),
            format.format(height),
            ha = "center",
            va = "bottom",
            color = textColor,
            fontsize = 14
        )

    return bars


def plot(contents: list[Result],
         quantityNames: list[str] = None,
         labelSizes: tuple[int,int] = (16, 16),
         labelColors: tuple[list[float],list[float]] = (colors["black"], colors["black"])) -> None:

    numberOfSubplots = len(contents)
    grid = (int(math.sqrt(numberOfSubplots)), 1)
    if grid[0]:
        grid = (grid[0], math.ceil(numberOfSubplots / grid[0]))

    figure, axisSets = pyplot.subplots(grid[0], grid[1], num = "Benchmark Result")
    if 1 < len(contents):
        axisSets = [axisSets]
    axisSets = numpy.ravel(axisSets)

    colorKeys = list(iter(colors.keys()))
    for results, axes in zip(contents, axisSets):
        labels = [case.name for case in results.cases]
        positions = numpy.arange(len(labels))
        for quantityIndex, quantityName in enumerate(quantityNames):
            colorKey = colorKeys[quantityIndex % len(colorKeys)]
            color = colors[colorKey]
            textColor = colors["gray"] if colorKey == "black" else colors["black"]
            graph = plotQuantity(
                results,
                quantityName,
                axes,
                barCount = len(quantityNames),
                barIndex = quantityIndex,
                positions = positions,
                color = color,
                textColor = textColor)
        axes.set_xticks(positions)
        axes.set_xticklabels(labels, rotation = 45, ha = "right")

        # Set y axis scale
        yRange = [float("inf"), float("-inf")]
        for rectangle in graph:
            height = rectangle.get_height()
            if height < yRange[0]:
                yRange[0] = height
            if yRange[1] < height:
                yRange[1] = height
        if yRange[0] and yRange[1] and 10 < yRange[1] / yRange[0]:
            axes.set_yscale("log")

    axes.grid(visible = True, axis = "y")

    figure.tight_layout()
    pyplot.show()
