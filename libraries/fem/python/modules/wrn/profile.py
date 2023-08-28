def run(*args, **kwargs):

    # --- Internal Imports
    from .heterogeneous_waveform_relaxation import run as target

    # --- STD Imports ---
    import pathlib
    import json
    import itertools
    import contextlib
    import io
    import re

    # Get default parameters
    with open(pathlib.Path(__file__).absolute().parent / "parameters.json", "r") as file:
        parameters = json.load(file)

    # Define parameter space
    parameterSpace = {
        "baseStep" : [parameters["baseStep"], parameters["baseStep"] / 10.0],
        "stepFrequencies" : [[1<<i_exponent, 1] for i_exponent in range(8)],
        "numberOfWindows" : [1, 5]
    }

    # Clear logs and init
    outputFilePath = pathlib.Path("wrn_profile.csv")
    with open(outputFilePath, "w") as file:
        file.write(",".join(parameterSpace.keys()))
        file.write("\n")
    print(*tuple(parameterSpace.keys()))

    outputRegex = re.compile("Elapsed time: (.*)")
    outputParser = lambda output: outputRegex.findall(output)[0]

    args = ["profile", "--quiet"]

    # Loop thorugh configurations
    for configuration in itertools.product(*list(parameterSpace.values())):
        # Apply configuration
        for key, value in zip(parameterSpace.keys(), configuration):
            parameters[key] = value

        # Run target with current configuration
        stream = io.StringIO()
        with contextlib.redirect_stdout(stream):
            target(args, parameters = parameters)

        # Write to output file
        elapsed = float(outputParser(stream.getvalue()))
        print(configuration, elapsed)
        with open(outputFilePath, "a") as file:
            file.write(",".join([str(value) for value in configuration]))
            file.write(f",{elapsed}")
            file.write("\n")


if __name__ == "__main__":
    run()
