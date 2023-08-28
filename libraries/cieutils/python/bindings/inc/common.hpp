#ifndef CIE_UTILS_PYTHON_BINDINGS_COMMON
#define CIE_UTILS_PYTHON_BINDINGS_COMMON

// --- External Includes ---
#include "packages/exceptions/inc/exception.hpp"
#include "pybind11/pybind11.h"

// --- Utility Includes ---
#include "packages/types/inc/types.hpp"
#include "packages/macros/inc/exceptions.hpp"


namespace cie {


template <class T>
void addBindings(Ref<pybind11::module_> r_module);


void addCMakeVariables(Ref<pybind11::module_> r_module);


} // namespace cie


#endif
