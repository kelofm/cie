"""Objects responsible for converting the output of benchmark tools into a common format."""

# --- Internal Imports ---
from .Result import Result

# --- SLD Imports ---
import abc
import io
import re
import json
import pathlib


class Parser(metaclass = abc.ABCMeta):

    @abc.abstractproperty
    def name(self) -> str:
        return ""

    @abc.abstractmethod
    def __call__(self, stream: io.TextIOBase) -> list[Result]:
        pass


class InternalJSONParser(Parser):
    """Parses a stream generated from serializing a list of @ref Result objects."""

    @property
    def name(self) -> str:
        return "json"

    def __call__(self, stream: io.TextIOBase) -> list[Result]:
        contents = json.load(stream)

        # The contents should either be a list of dicts or a single dict
        if isinstance(contents, dict):
            return [Result.deserialize(contents)]
        elif isinstance(contents, list):
            return [Result.deserialize(map) for map in contents]
        else:
            raise RuntimeError(f"Expecting a map or list of maps, but got {type(contents)}")


class GoogleTextParser(Parser):
    """Parses google benchmark text output and packs them into a single @ref Result object."""

    @property
    def name(self) -> str:
        return "google-text"

    def __init__(self):
        self.__regex = self.__makeRegex()

    @staticmethod
    def __makeRegex() -> re.Pattern:
        floatPattern = R"-?(?:(?:(?:[1-9][0-9]*)(?:\.[0-9]*)?)|(?:0(?:\.[0-9]*)?))(?:[eE][\+-]?[0-9]+)?"
        name = R".*"
        time = floatPattern
        timeUnit = R"\w+[^ ]"
        cpu = floatPattern
        cpuUnit = R"\w+[^ ]"
        iterationCount = R"[1-9][0-9]*"
        memoryThroughput = floatPattern
        memoryThroughputUnit = R"(?:k|M|G|T)?B\/s"

        pattern = f"^({name}) +({time}) ({timeUnit}) +({cpu}) ({cpuUnit}) +({iterationCount})"
        pattern += f" *(?:({memoryThroughput})({memoryThroughputUnit}))?.*$"
        return re.compile(pattern)

    def __call__(self, stream: io.TextIOBase) -> list[Result]:
        Case = Result.Case
        Quantity = Case.Quantity

        results = Result()
        for line in stream.readlines():
            line = line.strip()
            match = self.__regex.match(line)
            if match:
                results.addCase(Case(
                    name = match.group(1).strip(),
                    time = Quantity(float(match.group(2)), match.group(3)),
                    cpu = Quantity(float(match.group(4)), match.group(5)),
                    iterations = Quantity(float(match.group(6)), "")
                ))

        return [results]


class CIEProfilerParser(Parser):

    @property
    def name(self) -> str:
        return "cie-profiler"

    def __call__(self, stream: io.TextIOBase) -> list[Result]:
        Case = Result.Case
        Quantity = Case.Quantity

        root: dict = json.load(stream)
        if not all (key in root for key in ("results", "meta")):
            raise ValueError(f"Expecting JSON with 'results' and 'meta' but got one with the following keys\n: {root.keys()}")

        # Here's what the imported structure should look like:
        # {
        #   "meta" : {
        #       "name" : "total",
        #       "time" : some_integer,
        #       "timeUnit" : "time_unit_as_string",
        #       ...
        #   },
        #   "results" : [
        #       {
        #           "callCount" : some_integer,
        #           "file"      : "file_path_as_string",
        #           "function"  : "full_funciton_name_as_string",
        #           "line"      : line_number_as_integer,
        #           "time"      : some_integer
        #       },
        #       ...
        #   ]
        # }
        #
        # Case Name:
        # Parsing the function name is pretty much impossible without
        # parsing C++ because the entire function signature is in there
        # with all its potential template and namespace glory, so a
        # combination of the file path's stem and the line number will
        # make do as case name.
        contents: list[dict] = root["results"]
        meta = root["meta"]
        timeUnit = meta["timeUnit"].replace("µ", "u") # <== C++ prints µs

        # Add the total time to the contents
        results = Result()
        results.addCase(Case(
            name = meta["name"],
            time = Quantity(float(meta["time"]), timeUnit)
        ))

        # Loop through the contents and add them to the output
        for item in contents:
            results.addCase(Case(
                name = str(pathlib.Path(item["file"]).stem) + (f"_{item['line']}" if "line" in item else ""),
                time = Quantity(float(item["time"]), timeUnit)
            ))

        return [results]
