#pragma once

// --- Internal Includes ---
#include "poisson2D/definitions.hpp"

// --- FEM Includes ---
#include "packages/graph/inc/Assembler.hpp"

// --- Utility Includes ---
#include "packages/concurrency/inc/ThreadPoolBase.hpp"
#include "packages/commandline/inc/ArgParse.hpp"


namespace cie::fem {


void solve(
    CSRWrapper lhs,
    std::span<Scalar> solution,
    std::span<Scalar> rhs,
    Ref<const Assembler> rAssembler,
    Ref<mp::ThreadPoolBase> rThreads,
    Ref<const utils::ArgParse::Results> rArguments);


} // namespace cie::fem
