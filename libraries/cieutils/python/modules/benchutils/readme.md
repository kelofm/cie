# Description

```benchutils``` is a simple python module for parsing, plotting, and archiving the results of benchmarks / profiling runs.

# Interface

```benchutils``` can be executed as a module, instructed to read and parse a file or ```stdin```, filter the parsed data, plot it, and write it to ```stdout``` or a file. Arguments are parsed by ```argparse``` so you can take a peek at all options, but here's an overview:

- **```-i --input FILE```** : specify an input file to read, or ```-``` to read from ```stdin```
- **```-p --parser PARSER```** : name of the parser to use on the input. Currently, ```benchutils```' terrible internal format (```json```) and the text output of [google benchmark](https://github.com/google/benchmark) (```google-text```) are implemented. The default is ```google-text```.
- **```-o --output FILE```** : specify an output file to write to, or ```-``` for ```stdout```.
- **```--overwrite```** : overwrite the output file if it already exists. The default behaviour is to append the results.
- **```--plot NAME [NAMES ...]```** : plot a measured quantity of each case that has it. Plotting multiple quantities is *possible*.
- **```-e --external-modules MODULE [MODULES ...]```** : name of an external module to include build information from. The module must have an ```info``` variable storing its data as ```dict[str,str]```. The module must be on ```PYTHONPATH```.
- **```-f --filter ARGUMENT [ARGUMENTS ...]```** : apply a filter on the parsed results to reduce the amount of data in the working set. Filter arguments are generated dynamically, based on what was parsed from the input. In any case, each argument must have exactly one string argument, representing a regex pattern. Cases that match all provided regexes are **kept** in the working set and all others are discarded (yes, I know it's exactly the inverse of what a filter does 😢).

Examples for filter arguments:

- **```--filter="--OS-system Linux --cases-name parallel"```** keeps all results from cases that were run on a Linux system and have ```parallel``` in their names.
- **```-f "--python-compiler \"GCC (12)|(11)\""```** keeps results from cases that were run on a system that has its python interpreter compiled with GCC 12 or 11.

Again, filter arguments are parsed using ```argparse```, so passing a ```-h``` will help you out.

Example for a full command:
```
python -m benchutils --input some/path/some_archive.json --parser json --filter "--OS-release MANJARO" --plot time cpu
```

The command above reads an existing file at ```some/path/some_archive.json``` and parse it assuming it was written in the internal ```json``` format. Then all cases are selected that were run on [Manjaro](https://manjaro.org/), then ```time``` and ```cpu``` are plotted from cases that have them.

# ToDo

- Make plots nicer :(
- Implement some logic that groups runs / cases together, finds out what the most meaningful differences between them are, and plots + labels them accordingly.
- Implement some other parsers (eg.: google benchmark's json and csv output).
