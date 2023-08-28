# --- STD Imports ---
import time
import contextlib
import argparse
import pathlib


@contextlib.contextmanager
def TimedContext() -> float:
    begin = time.perf_counter()
    yield lambda: time.perf_counter() - begin


def parseArguments(argv: list) -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser("wrn")
    parser.add_argument("-i",
                        "--input-path",
                        dest = "inputPath",
                        type = pathlib.Path,
                        default = pathlib.Path(__file__).absolute().parent / "parameters.json")
    parser.add_argument("-q",
                        "--quiet",
                        dest = "quiet",
                        action = "store_const",
                        default = False,
                        const = True,
                        help = "Suppress status reports")
    parser.add_argument("-p",
                        "--plot",
                        dest = "plot",
                        action = "store_const",
                        default = False,
                        const = True,
                        help = "plot selected DoFs")
    parser.add_argument("-w",
                        "--write",
                        dest = "write",
                        action = "store_const",
                        default = False,
                        const = True,
                        help = "write frames")
    return parser.parse_args(argv[1:])
