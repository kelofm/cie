// --- External Includes ---
#include "pybind11/pybind11.h"

// --- Utility Includes ---
#include "packages/types/inc/types.hpp"


namespace cie::fem::maths {


void makeFEMMathsBindings(Ref<pybind11::module_> rModule);


} // namespace cie::fem::maths
