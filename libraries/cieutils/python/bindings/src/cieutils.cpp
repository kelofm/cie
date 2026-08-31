// --- External Includes ---
#include "pybind11/pybind11.h"

// --- Utility Includes ---
#include "version.hpp"
#include "python/bindings/inc/common.hpp"


PYBIND11_MODULE(cieutils_python_bindings, rModule) {
    rModule.doc() = "CiE Utilities";
    auto propertyModule = rModule.def_submodule("properties");
    cie::addBindings<cie::Version>(propertyModule);
    cie::addCMakeVariables(propertyModule);
} // PYBIND_MODULE
