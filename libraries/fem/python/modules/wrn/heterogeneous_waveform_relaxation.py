def run(args, parameters: dict = {}):
    # --- External Imports ---
    import numpy

    # --- Internal Imports ---
    from .time_integration import NewmarkTimeIntegrator
    from .model import Model, Solution
    from .relax import HeterogeneousJacobiRelaxation
    from .windowing import WindowSequence
    from .utilities import TimedContext, parseArguments

    # --- STD Imports ---
    import json


    # CLI
    args = parseArguments(args)


    # Parse system parameters
    if not parameters:
        with open(args.inputPath, "r") as file:
            parameters = json.load(file)


    masses          = parameters["masses"]
    dampings        = parameters["dampings"]
    stiffnesses     = parameters["stiffnesses"]
    baseStep        = parameters["baseStep"]
    stepFrequencies = [int(item) for item in parameters["stepFrequencies"]]
    timeEnd         = parameters["timeEnd"]
    numberOfWindows = parameters["numberOfWindows"]
    watchDoFs       = parameters["watchDoFs"]


    model = Model(numpy.array([[masses[0], 0.0], [0.0, masses[1]]]),                                        # mass matrix
                  numpy.array([[dampings[0], 0.0], [0.0, dampings[1]]]),                                    # damping matrix
                  numpy.array([[sum(stiffnesses), -stiffnesses[1]], [-stiffnesses[1], stiffnesses[1]]]),    # stiffness matrix
                  lambda step, time: numpy.array([0.0, 0.0 if time < timeEnd / 2.0 else -1e1]),             # load vector
                  [[0.0, 1e0], [0.0, 0.0]])                                                                 # initial conditions [i_displacement, i_velocity]

    # ------------------------------------------------------------------------ #

    # Computed values
    stepSizes = [baseStep * frequency for frequency in stepFrequencies]
    numberOfSteps = [round(timeEnd / dt) + 1 for dt in stepSizes]
    solution = Solution(numpy.linspace(0.0, timeEnd, max(numberOfSteps)), len(masses), stepFrequencies)


    # Traditional solve on the fast scale
    if args.plot:
        referenceSolution = solution.copy()
        referenceTimeIntegrator = NewmarkTimeIntegrator(min(stepSizes), len(solution.timeSamples))
        referenceDisplacements, referenceVelocities, referenceAccelerations = referenceTimeIntegrator.integrate(model)
        for i_variable in range(solution.numberOfVariables):
            referenceSolution.set(referenceDisplacements[:,i_variable], i_variable, dense = True)
            referenceSolution.set(referenceVelocities[:,i_variable], i_variable, derivative = 1, dense = True)
            referenceSolution.set(referenceAccelerations[:,i_variable], i_variable, derivative = 2, dense = True)


    # Graphics setup
    class SaveFig:
        def __init__(self):
            self.__counter = 0
        def __call__(self):
            pyplot.savefig(f"wrn_{self.__counter:02}.png",
                           dpi = 800)
            self.__counter += 1
    saveFig = SaveFig()

    if args.plot or args.write:
        from matplotlib import pyplot
        figure, axes = pyplot.subplots()
        referenceLines = []
        relaxedLines = []
        legend = []
        colors = [
            [component / 255.0 for component in (48, 112, 179)], # <== TUM web blue
            [component / 255.0 for component in (227, 114, 34)], # <== TUM orange
            [component / 255.0 for component in (162, 173, 0)]   # <== TUM green
        ]
        for i_variable in watchDoFs:
            for line in axes.plot(solution.timeSamples,
                                  referenceSolution.get(i_variable, dense = True),
                                  alpha = 0.5,
                                  color = colors[i_variable % len(colors)],
                                  linestyle = "dashed"):
                referenceLines.append(line)
                legend.append(f"Ref DoF #{i_variable}")
            for line in axes.plot(solution.getTimeSamples(i_variable),
                                  solution.get(i_variable),
                                  color = colors[i_variable % len(colors)]):
                relaxedLines.append(line)
                legend.append(f"DoF #{i_variable}")
        axes.legend(legend)
        axes.xaxis.label.set_text("t")
        axes.yaxis.label.set_text("u")

        if args.plot:
            pyplot.show(block = False)

        if args.write:
            saveFig()
    else:
        referenceLines = None
        relaxedLines = None


    windows = WindowSequence(model, solution, numberOfWindows)


    def hook(i_window: int, window: WindowSequence.Window, i_relax: int = 0, norm: float = 0.0) -> None:
        if not args.quiet:
            print(f"Window {i_window} | Relaxation {i_relax} | Residual: {norm:.2E}")

        if args.plot or args.write:
            windows.updateSolution()
            for i_variable, line in zip(watchDoFs, relaxedLines):
                line.set_ydata(solution.get(i_variable))
            figure.canvas.draw()
            figure.canvas.flush_events()
            if args.write:
                saveFig()


    with TimedContext() as getElapsed:
        for i_window, window in enumerate(windows):
            relaxation = HeterogeneousJacobiRelaxation(window.model, window.solution)
            relaxation.solve(hook = lambda i_relax, norm: hook(i_window, window, i_relax, norm))
        print(f"Elapsed time: {getElapsed():.3F}")


    if args.plot:
        pyplot.show()


if __name__ == "__main__":
    import sys
    run(sys.argv)
