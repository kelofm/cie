#pragma once

// --- FEM Includes ---
#include "packages/io/inc/VTKHDF.hpp"


namespace cie::io {


template <fem::DiscretizationLike TMesh>
void VTKHDF::Output::operator()(Ref<const TMesh>) {}


} // namespace cie::io
