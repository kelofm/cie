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


def findMinPrefix(cases: list[Result.Case], quantityName: str):
    minPrefix = (None, None)
    prefixes = set()

    if isinstance(cases, list):
        for case in cases:
            prefixes.add(case[quantityName].getPrefix())
    elif isinstance(cases, Result.Case):
        for case in cases:
            prefixes.add(case[quantityName].getPrefix())
    else:
        return ValueError(f"Unexpected input type: {type(cases)}")

    for prefix in prefixes:
        exponent = Result.Case.Quantity.getPrefixMap()[prefix]
        if minPrefix[1] == None or minPrefix[2] < exponent:
            minPrefix = (prefix, exponent)

    return minPrefix[0]


def homogenizeUnits(cases: list[Result.Case], quantityName: str, commonPrefix: str = "") -> None:
    if not commonPrefix:
        commonPrefix = findMinPrefix(cases, quantityName)

    for case in cases:
        case[quantityName].setPrefix(commonPrefix)


def addBars(cases: list[Result.Case],
            quantityName: str,
            axes: pyplot.Axes,
            color: list[float] = colors["TUMBlue"],
            prefix: str = None,
            label: str = "",
            positions: numpy.array = None,
            barCount: int = 1,
            barIndex: int = 0,
            totalWidth: float = 0.8,
            format = R"{:.2E}") -> typing.Any:
    if not prefix:
        prefix = findMinPrefix(cases, quantityName)
    homogenizeUnits(cases, quantityName, commonPrefix = prefix)

    if not label:
        label = quantityName

    # Collect values
    values = [case[quantityName].getValue() for case in cases]

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
        textColor = [0.0, 0.0, 0.0, 1.0] if barIndex % 2 else [1.0, 1.0, 1.0, 1.0]
        axes.text(position + width / 2.0,
                  height * ((barIndex + 1) / barCount),
                  format.format(height),
                  ha = "center",
                  va = "bottom",
                  color = textColor,
                  fontsize = 14)

    axes.set_title(f"{quantityName} [{cases[0][quantityName].unit}]")
    return bars


def plot(contents: list[Result],
         quantityNames: list[str] = [],
         labelSizes: tuple[int,int] = (16, 16),
         labelColors: tuple[list[float],list[float]] = (colors["black"], colors["black"])) -> None:
    # Find and group matching cases in the contents by name
    groupedCases: dict[str,list[Result.Case]] = {}
    caseNames: set[str] = [] # <== not built using a set to preserve original order
    for results in contents:
        for case in results.cases:
            if not case.name in caseNames:
                caseNames.append(case.name)

    for caseName in caseNames:
        for results in contents:
            case = next((c for c in results.cases if c.name == caseName), None)
            groupedCases.setdefault(caseName, []).append(case)

    numberOfSubplots = len(quantityNames)
    grid = (int(math.sqrt(numberOfSubplots)), 1)
    if grid[0]:
        grid = (grid[0], math.ceil(numberOfSubplots / grid[0]))

    figure, axisSets = pyplot.subplots(grid[0], grid[1], num = "Benchmark Results")
    if 1 < len(contents):
        axisSets = [axisSets]
    axisSets = numpy.ravel(axisSets)

    colorKeys = list(iter(colors.keys()))
    labels = caseNames
    positions = numpy.arange(len(caseNames))
    for quantityName, axes in zip(quantityNames, axisSets):
        for i_result in range(len(contents)):
            colorKey = colorKeys[i_result % len(colorKeys)]
            color = colors[colorKey]
            textColor = colors["gray"] if colorKey == "black" else colors["black"]
            graph = addBars([cases[i_result] for cases in groupedCases.values()],
                            quantityName,
                            axes,
                            barCount = len(contents),
                            barIndex = i_result,
                            positions = positions,
                            color = color)
        axes.set_xticks(positions)
        axes.set_xticklabels(caseNames, rotation = 45, ha = "right")

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
