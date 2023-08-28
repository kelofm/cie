// --- External Includes ---
#include "pybind11/pybind11.h"

// --- FEM Includes ---
#include "python/bindings/maths.hpp"


// ---------------------------------------------------------
// MODULE DEFINITION
// ---------------------------------------------------------
PYBIND11_MODULE(fem_python_bindings, pybindModule)
{
    cie::fem::maths::makeFEMMathsBindings(pybindModule);
} // PYBIND11_MODULE __fem
