def loadBindings(
    libraryName: str,
    namespace: dict) -> None:
        # --- STD Imports ---
        import pathlib as __pathlib
        from importlib import util as __importutil
        import sys as __sys
        import re as __re
        from typing import Optional as __Optional

        # Find binary package
        __scriptDir = __pathlib.Path(__file__).absolute().parent
        __pythonBindingsPattern: __re.Pattern = __re.compile(f"({libraryName}_python_bindings).*")
        __libDir = __scriptDir.parent
        __pythonBindings = __libDir.glob("*")

        if not __pythonBindings:
            raise FileNotFoundError(f"Could not find bindings for cieutils in ${__libDir}")

        for __maybeBindings in __pythonBindings:
            __match: __Optional[__re.Match] = __pythonBindingsPattern.match(__maybeBindings.name)
            if __match:
                # Import the binary module
                __spec = __importutil.spec_from_file_location(__match.group(1), str(__maybeBindings))
                __module = __importutil.module_from_spec(__spec)
                __sys.modules[__spec.name] = __module
                __spec.loader.exec_module(__module)

                # Import everything from the binary module
                # and fill the current namespace
                for __attribute in dir(__module):
                    if not __attribute.startswith("_"):
                        namespace[__attribute] = getattr(__module, __attribute)


loadBindings("cieutils", globals())
