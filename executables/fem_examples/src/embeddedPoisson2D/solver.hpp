#pragma once

// --- Internal Includes ---
#include "embeddedPoisson2D/definitions.hpp"

// --- FEM Includes ---
#include "packages/graph/inc/Assembler.hpp"

// --- Utility Includes ---
#include "packages/concurrency/inc/ThreadPoolBase.hpp"
#include "packages/io/inc/json.hpp"


namespace cie::fem {


void solve(
    linalg::CSRView<Scalar,int> lhs,
    std::span<Scalar> solution,
    std::span<Scalar> rhs,
    Ref<Assembler> rAssembler,
    Ref<mp::ThreadPoolBase> rThreads,
    Ref<const cie::io::JSONObject> rConfiguration);


} // namespace cie::fem
