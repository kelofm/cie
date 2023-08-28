def run(args: list, parameters: dict = {}):
    # --- External Imports ---
    import numpy

    # --- Internal Imports ---
    import cie.wrn as wrn
    from cie.wrn.time_integration import NewmarkTimeIntegrator
    from cie.wrn.model import Model, Solution
    from cie.wrn.relax import HomogeneousRelaxation
    from cie.wrn.windowing import WindowSequence
    from cie.wrn.utilities import TimedContext, parseArguments

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
    stepFrequency   = min(parameters["stepFrequencies"])
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
    stepSize = baseStep * stepFrequency
    numberOfSteps = round(timeEnd / stepSize) + 1
    solution = Solution(numpy.linspace(0.0, timeEnd, numberOfSteps), len(masses), [stepFrequency for dummy in masses])

    # Traditional solve on the fast scale
    referenceSolution = solution.copy()
    referenceTimeIntegrator = NewmarkTimeIntegrator(stepSize, numberOfSteps)
    referenceDisplacements, referenceVelocities, referenceAccelerations = referenceTimeIntegrator.integrate(model)
    for i_variable in range(solution.numberOfVariables):
        referenceSolution.set(referenceDisplacements[:,i_variable], i_variable, dense = True)
        referenceSolution.set(referenceVelocities[:,i_variable], i_variable, derivative = 1, dense = True)
        referenceSolution.set(referenceAccelerations[:,i_variable], i_variable, derivative = 2, dense = True)


    # Graphics setup
    if args.plot:
        from matplotlib import pyplot
        figure, axes = pyplot.subplots()
        referenceLines = []
        relaxedLines = []
        legend = []
        for i_variable in watchDoFs:
            for line in axes.plot(solution.getTimeSamples(i_variable), referenceSolution.get(i_variable)):
                referenceLines.append(line)
                legend.append(f"Ref DoF #{i_variable}")
            for line in axes.plot(solution.getTimeSamples(i_variable), solution.get(i_variable)):
                relaxedLines.append(line)
                legend.append(f"DoF #{i_variable}")
        axes.legend(legend)
        pyplot.show(block = False)
    else:
        referenceLines = None
        relaxedLines = None


    windows = WindowSequence(model, solution, numberOfWindows)


    def hook(i_window: int, window: WindowSequence.Window, i_relax: int = 0, norm: float = 0.0) -> None:
        print(f"Window {i_window} | Relaxation {i_relax} | Residual: {norm:.2E}")
        windows.updateSolution()

        if args.plot:
            for i_variable, line in enumerate(relaxedLines):
                line.set_ydata(solution.get(i_variable))
            figure.canvas.draw()
            figure.canvas.flush_events()


    with TimedContext() as getElapsed:
        for i_window, window in enumerate(windows):
            relaxation = HomogeneousRelaxation(window.model, window.solution)
            relaxation.solve(hook = lambda i_relax, norm: hook(i_window, window, i_relax, norm))
        print(f"Elapsed time: {getElapsed():.3F}")

    if args.plot:
        pyplot.show()


if __name__ == "__main__":
    import sys
    run(sys.argv)
