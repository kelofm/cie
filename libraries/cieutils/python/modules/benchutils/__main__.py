# --- Internal Includes ---
from .Result import Result
from .cli import getParsers, makeArgumentParser, makeFilterArgumentParser, filterResult

# --- STD Imports ---
import pathlib
import sys
import json
import importlib


# Parse primary arguments
argumentParser = makeArgumentParser()
arguments = argumentParser.parse_args(sys.argv[1:])

# Parse input
benchmarkParser = getParsers()[arguments.parser]
if arguments.input == "-":
    contents = benchmarkParser(sys.stdin)
else:
    with open(pathlib.Path(arguments.input), 'r') as file:
        contents = benchmarkParser(file)

for moduleName in arguments.externalModules:
    module = importlib.import_module(moduleName)
    infoDict = module.info

    for results in contents:
        results.addExternalInfo(moduleName, infoDict)

# Filter cases
if arguments.filter:
    filterArguments = makeFilterArgumentParser(contents).parse_args(
        argument for argument in arguments.filter.split(' ') if argument
    )
    filterResult(contents, filterArguments)

# Write output
if arguments.output:
    output = arguments.output
    outputFile = pathlib.Path(output)

    if outputFile.is_file() and not arguments.overwrite:
        with open(outputFile, 'r') as file:
            archive = json.load(file)
    else:
        archive = []

    archive += [results.serialize() for results in contents]

    if output == '-':
        for results in archive:
            print(json.dumps(results, indent = 4))
    elif output:
        with open(outputFile, 'w') as file:
            json.dump(archive, file, indent = 2)
    del archive

# Plot
if arguments.plot:
    from .graphics import plot
    plot(contents, arguments.plot)
