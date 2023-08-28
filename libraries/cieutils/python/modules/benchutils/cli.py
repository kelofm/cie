# --- Internal Imports ---
from .Result import Result
from .Parser import Parser

# --- STD Imports ---
import argparse
import re


# Constrained argument options
def getParsers() -> dict[str,Parser]:
    return dict([(cls().name, cls()) for cls in Parser.__subclasses__()])


def makeArgumentParser() -> argparse.ArgumentParser:
    argumentParser = argparse.ArgumentParser("benchutils")
    benchmarkParsers = getParsers()

    argumentParser.add_argument("-i",
                                "--input",
                                dest = "input",
                                type = str,
                                default = "-",
                                help = "Specify the input to be parsed. Can be a path to a file, or '-' for stdin")

    argumentParser.add_argument("-p",
                                "--parser",
                                dest = "parser",
                                type = str,
                                choices = list(benchmarkParsers.keys()),
                                default = next(iter(benchmarkParsers.keys())))

    argumentParser.add_argument("-o",
                                "--output",
                                dest = "output",
                                nargs = '?',
                                default = '-',
                                help = "Specify the output path, or '-' for stdout (default)")

    argumentParser.add_argument("--overwrite",
                                dest = "overwrite",
                                action = "store_const",
                                const = True,
                                default = False,
                                help = "Overwrite the output file")

    argumentParser.add_argument("-f",
                                "--filter",
                                dest = "filter",
                                type = str,
                                default = "",
                                help = "Apply regex-based filters to the parsed results. Options are generated dynamically => pass -h as an argument to see what's available.")

    argumentParser.add_argument("--plot",
                                dest = "plot",
                                nargs = '+',
                                default = [])

    argumentParser.add_argument("-e",
                                "--external-modules",
                                dest = "externalModules",
                                nargs = '+',
                                default = [],
                                help = "Import external modules and include their 'info' dicts")

    return argumentParser


def makeFilterArgumentParser(contents: list[Result]) -> argparse.ArgumentParser:
    # Collect valid arguments
    regexTargets = set()
    for results in contents:
        # Everything except 'cases' is expected to be a string-string map
        for key in results.keys():
            if key != "cases":
                for subKey in results[key].keys():
                    regexTargets.add(f"{key}-{subKey}")

        # 'cases' has special treatment because it's a list of Case objects
        for case in results.cases:
            for subKey in case.keys():
                regexTargets.add(f"cases-{subKey}")

    # Build filter parser
    filterParser = argparse.ArgumentParser("benchutils-filter")
    for regexTarget in regexTargets:
        filterParser.add_argument(f"--{regexTarget}",
                                  dest = f"{regexTarget.replace('-', '_')}",
                                  type = lambda pattern: re.compile(pattern))
    return filterParser


def filterResult(contents: list[Result], filterArguments: argparse.Namespace) -> None:
    # Parse and apply filters
    infoPredicates = []
    casePredicates = []
    for regexTarget, regex in vars(filterArguments).items():
        if regex != None:
            objectPath = regexTarget.split('_')
            if objectPath[0] == "cases":
                casePredicates.append(lambda case, objectPath = objectPath, regex = regex: not regex.match(getattr(case, objectPath[1])))
            else:
                infoPredicates.append(lambda results, objectPath = objectPath, regex = regex: regex.match(
                    getattr(results, objectPath[0])[objectPath[1]]
                ))

    # Apply info predicates
    contents[:] = [results for results in contents if all(predicate(results) for predicate in infoPredicates)]

    # Apply case predicates and delete empty results
    for results in contents:
        for predicate in casePredicates:
            results.filterCases(predicate)
    contents[:] = [results for results in contents if results.cases]
