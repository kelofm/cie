"""A common format used throughout this module."""

# --- STD Imports ---
import platform
import typing
import functools
import re
import subprocess
import os


def tryFetchFromExternal(callable: typing.Callable[[], str]) -> str:
    @functools.wraps(callable)
    def wrapped() -> str:
        try:
            return callable()
        except ModuleNotFoundError:
            return "unavailable"
    return wrapped


class SystemInfo:
    """ @brief A dict restricted to "OS", "hardware", and "python"."""

    def __init__(self, inputMap: dict = {}):
        """ @brief Reconstruct a SystemInfo instance from a dict or initialize one with the current system's info.
            @details Only "OS", "hardware", and "python" are copied from @a inputDict and the rest of the entries
                     are ignored. Missing entries are added using the current system's info. The internal dict has
                     the following structure: @code {
                        "OS" : {
                            "system" : "",
                            "release" : "",
                            "version" : ""
                        },
                        "hardware" : {
                            "machine" : "",
                            "processor" : "",
                            "memory" : ""
                        },
                        "python" : {
                            "build" : "",
                            "compiler" : "",
                            "version" : ""
                        }
                     } @endcode
            @note Fetching the size of the system's main memory requires the 'psutil' package.
            @note Properties that could not be fetched are replaced by "unavailable".
        """
        self.__map = {key : value for key, value in inputMap.items() if key in ("OS", "hardware", "python")} if inputMap else {}
        self.__map.setdefault("OS", self.__getOS())
        self.__map.setdefault("hardware", self.__getHardware())
        self.__map.setdefault("python", self.__getPython())

    def serialize(self) -> dict:
        return self.__map.copy()

    @staticmethod
    def deserialize(map: dict) -> "SystemInfo":
        return SystemInfo(inputMap = map)

    def keys(self) -> list[str]:
        return list(iter(self.__map.keys()))

    def __contains__(self, item: str) -> bool:
        return item in self.__map

    def __getitem__(self, key: str) -> dict:
        return self.__map[key]

    @staticmethod
    def __getOS() -> dict[str,str]:
        OS = platform.system_alias(platform.system(), platform.release(), platform.version())
        return {
            "system" : OS[0],
            "release" : OS[1],
            "version" : OS[2]
        }

    @staticmethod
    def __getHardware() -> dict[str,str]:
        return {
            "machine" : platform.machine(),
            "processor" : SystemInfo.__getCPU(),
            "memory" : SystemInfo.__getMemory()
        }

    @staticmethod
    def __getPython() -> dict:
        return {
            "build" : " ".join(platform.python_build()),
            "compiler" : platform.python_compiler(),
            "version" : ".".join(platform.python_version_tuple())
        }

    @tryFetchFromExternal
    @staticmethod
    def __getMemory() -> str:
        import psutil
        return str(f"{psutil.virtual_memory().total / (1<<30):.3}") + " GB"

    @staticmethod
    def __getCPU() -> str:
        if platform.system() == "Linux":
            wallOfText = subprocess.check_output("cat /proc/cpuinfo", shell = True).decode().strip()
            for line in wallOfText.split('\n'):
                if "model name" in line:
                    return re.sub(".*model name.*:", "", line, 1)

        elif platform.system() == "Windows":
            return platform.processor()

        elif platform.system() == "Darwin":
            os.environ['PATH'] = os.environ['PATH'] + os.pathsep + '/usr/sbin'
            command ="sysctl -n machdep.cpu.brand_string"
            return subprocess.check_output(command).strip()

        else:
            return "unavailable"


class Result:
    """ @brief Result is a @a dict -like object containing benchmark cases + system and package information.
        @details An instance of @ref Result is meant to represent all benchmark cases run on one system,
                 with one version of the test package and a specific set of dependencies. Package and
                 dependency versions can be passed as external information, which is empty by default.
    """

    class Case:
        """ @brief Class representing a benchmark case.
            @details @a Case behaves like a @a dict and always has a "name" key with an associated string.
                     The rest of the entries are parsed from keyword arguments passed to the constructor,
                     all of which must have Result.Case.Quantity as value.
        """

        class Quantity:
            """ @brief Class representing a value with an associated unit."""

            __prefixMap = {'T' : 12,
                           'G' : 9,
                           'M' : 6,
                           'k' : 3,
                           '': 0,
                           'm' : -3,
                           'u' : -6,
                           'n' : -9}

            __regex = re.compile(R"(-?(?:(?:(?:[1-9][0-9]*)(?:\.[0-9]*)?)|(?:0(?:\.[0-9]*)?))(?:[eE][\+-]?[0-9]+)?) ?([a-zA-Z_]+)?")

            def __init__(self, value: float, unit: str):
                self.__value = value
                self.__unit = unit

            def getValue(self, unit: str = None) -> float:
                if not unit:
                    return float(self.__value)
                else:
                    if self.isCompatibleWith(unit):
                        return self.__value * 10**(self.__parsePrefix(unit) - self.__parsePrefix(self.__unit))
                    else:
                        raise ValueError(f"{unit} is incompatible with {self.__unit}")

            def getCoefficient(self) -> float:
                return 10**(self.__prefixMap[self.getPrefix()])

            def getPrefix(self) -> str:
                return self.__getSingleUnitParts(self.__unit)[0]

            def setPrefix(self, prefix: str) -> None:
                if prefix not in self.__prefixMap:
                    raise ValueError(f"Invalid prefix: {prefix}")
                oldPrefix = self.__getSingleUnitParts(self.__unit)[0]
                self.__value *= 10**(self.__prefixMap[prefix] - self.__prefixMap[oldPrefix])
                self.__unit = prefix + self.__getSingleUnitParts(self.__unit)[1]

            def isCompatibleWith(self, unit: str) -> bool:
                return self.__getSingleUnitParts(self.__unit)[1] == self.__getSingleUnitParts(unit)[1]

            def serialize(self, unit: str = None) -> str:
                unit = unit if unit else self.__unit
                return f"{self.getValue(unit)} {unit if unit != '0' else ''}"

            @staticmethod
            def deserialize(string: str) -> "Result.Case.Quantity":
                match = Result.Case.Quantity.__regex.match(string)
                if match:
                    return Result.Case.Quantity(float(match.group(1)), match.group(2) if match.group(2) else "0")
                else:
                    raise ValueError(f"Failed to deserialize {string}")

            def __str__(self) -> str:
                return self.serialize()

            @property
            def unit(self) -> str:
                return self.__unit

            @staticmethod
            def getPrefixMap() -> dict[str,int]:
                return Result.Case.Quantity.__prefixMap

            @staticmethod
            def __parsePrefix(unit: str) -> int:
                """Parse the input unit and return the 10-based exponent of its SI coefficient"""
                unit = unit.split('/')
                if len(unit) == 1:
                    return Result.Case.Quantity.__parseSingleUnit(unit[0])
                elif len(unit) == 2:
                    return Result.Case.Quantity.__parseSingleUnit(unit[0]) - Result.Case.Quantity.__parseSingleUnit(unit[1])
                else:
                    raise ValueError(f"Unexpected unit: {unit}")

            @staticmethod
            def __parseSingleUnit(unit: str) -> int:
                # The unit may consist of 1 (base only) or more (prefix + base) characters
                return Result.Case.Quantity.__prefixMap[Result.Case.Quantity.__getSingleUnitParts(unit)[0]]

            @staticmethod
            def __getSingleUnitParts(unit: str) -> tuple[str,str]:
                """Return (prefix, base)"""
                if not unit:
                    raise ValueError("Empty unit")
                elif len(unit) == 1:
                    return "", unit
                else:
                    return unit[0], unit[1:]


        def __init__(self,
                     name: str,
                     **kwargs: dict[str,Quantity]):
            """ @brief Construct a @ref Case with a name and optional Result.Case.Quantity instances.
                @details All keyword arguments represent named quantities and are considered benchmark
                         output quantities (for example: 'time', 'cpuTime', or 'throughput').
            """
            self.__name = name
            self.__quantities = kwargs

            # Validate
            for key, value in self.__quantities.items():
                if not isinstance(value, Result.Case.Quantity):
                    raise ValueError(f"Expecting a Quantity but got {type(value)} for key {key} instead")

        @property
        def name(self) -> str:
            return self.__name

        def keys(self) -> list[str]:
            return ["name"] + list(iter(self.__quantities.keys()))

        def __getitem__(self, key: str) -> typing.Any:
            if key == "name":
                return self.__name
            elif key in self.__quantities:
                return self.__quantities[key]
            else:
                raise AttributeError(f"Case {self.name} does not have attribute {key}")

        def __getattr__(self, attributeName: str) -> typing.Any:
            return self[attributeName]

        def serialize(self, units: list[str] = None) -> dict:
            output = {"name" : self.name}
            for key, value in self.__quantities.items():
                output[key] = next((value.serialize(unit) for unit in units if value.isCompatibleWith(unit)), None if units else value.serialize())

            # Validate
            for key, value in output.items():
                if value == None:
                    raise ValueError(f"No unit in {units} is compatible with {value.unit} in {key}")

            return output

        @staticmethod
        def deserialize(map: dict) -> "Result.Case":
            return Result.Case(
                map["name"],
                **{key : Result.Case.Quantity.deserialize(value) for key, value in map.items() if key != "name"}
            )


    def __init__(self,
                 cases: list[Case] = None,
                 systemInfo: SystemInfo = None,
                 externalInfo: dict = None):
        """ @brief Construct benchmark results from the related cases, system info, and package + dependency info.
            @details Package and dependency information can be passed via the @a externalInfo input @a dict.
        """
        self.__cases = cases if cases else []
        self.__systemInfo = systemInfo if systemInfo else SystemInfo()
        self.__externalInfo = externalInfo if externalInfo else {}

    @property
    def info(self) -> SystemInfo:
        return self.__systemInfo

    @property
    def externalInfo(self) -> dict:
        return self.__externalInfo

    @property
    def cases(self) -> list[Case]:
        return self.__cases

    def addCase(self, case: Case) -> None:
        self.__cases.append(case)

    def addExternalInfo(self, key: str, value: dict) -> None:
        if not isinstance(key, str) or not isinstance(value, dict):
            raise ValueError(f"Invalid external type: {type(key)} : {type(value)}")
        self.__externalInfo[key] = value

    def filterCases(self, predicate: typing.Callable[["Result.Case"],bool]) -> None:
        self.__cases = [case for case in self.__cases if not predicate(case)]

    def keys(self) -> list[str]:
        return ["cases"] + list(iter(self.__systemInfo.serialize().keys())) + list(iter(self.__externalInfo.keys()))

    def hasKey(self, key: str) -> bool:
        return key in self.keys()

    def hasCase(self, caseName: str) -> bool:
        return next((True for case in self.__cases if case.name == caseName), False)

    def __getitem__(self, key: str) -> dict:
        if key == "cases":
            return self.__cases
        elif key in self.__systemInfo:
            return self.__systemInfo[key]
        elif key in self.__externalInfo:
            return self.__externalInfo[key]
        else:
            raise KeyError(f"Invalid key: {key}")

    def __getattr__(self, key: str) -> dict:
        return self[key]

    def serialize(self, units: list[str] = []) -> dict:
        output = {"cases" : [case.serialize(units) for case in self.__cases]}
        output.update(self.__systemInfo.serialize())
        output.update(self.__externalInfo)
        return output

    @staticmethod
    def deserialize(map: dict) -> "Result":
        """Input map must have 'results' key and be deserializable by @ref SystemInfo."""
        cases = [Result.Case.deserialize(case) for case in map["cases"]]
        systemInfo = SystemInfo.deserialize(map)
        externalInfo = {key : value for key, value in map.items() if key != "cases" and key not in systemInfo}
        return Result(cases = cases, systemInfo = systemInfo, externalInfo = externalInfo)
