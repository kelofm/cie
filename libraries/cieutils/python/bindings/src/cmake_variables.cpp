// --- External Includes ---
#include "pybind11/pybind11.h"
#include "pybind11/stl/filesystem.h"

// --- Utility Includes ---
#include "cmake_variables.hpp"
#include "python/bindings/inc/common.hpp"


namespace cie {


void addCMakeVariables(Ref<pybind11::module_> r_module)
{
    r_module.def("getInstallPath", &cie::getInstallPath)
            .def("getBinaryPath", &cie::getBinaryPath)
            .def("getSourcePath", &cie::getSourcePath)
            .def("getDataPath", &cie::getDataPath)
            .def("getOutputPath", &cie::getOutputPath)
            .def("getCompiler", &cie::getCompiler)
            .def("getCompilerVersion", &cie::getCompilerVersion)
            .def("getCompileFlags", &cie::getCompileFlags)
            .def("getCompilePlatform", &cie::getCompilePlatform)
            ;
}


} // namespace cie
