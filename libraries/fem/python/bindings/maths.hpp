// --- External Includes ---
#include "pybind11/pybind11.h"

// --- Utility Includes ---
#include "packages/types/inc/types.hpp"


namespace cie::fem::maths {


void makeFEMMathsBindings(Ref<pybind11::module_> r_module);


} // namespace cie::fem::maths
