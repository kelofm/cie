# This is the __init__ of the root cie module

# cieutils is guaranteed to be built
from . import utils as __cieutils

info = {
    "version" : str(__cieutils.properties.Version.local),
    "compiler" : f"{__cieutils.properties.getCompiler()} {__cieutils.properties.getCompilerVersion()}",
    "compilerFlags" : __cieutils.properties.getCompileFlags(),
    "compilePlatform" : __cieutils.properties.getCompilePlatform()
}
