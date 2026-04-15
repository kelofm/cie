#pragma once

// --- Internal Includes ---
#include "embeddedPoisson2D/CellData.hpp"
#include "embeddedPoisson2D/mesh.hpp"
#include "embeddedPoisson2D/integration.hpp"
#include "embeddedPoisson2D/constraints.hpp"

// --- Utility Includes ---
#include "packages/commandline/inc/ArgParse.hpp"


namespace cie::fem {


void postprocess(
    std::span<const Scalar> meshBase,
    std::span<const Scalar> meshLengths,
    linalg::CSRView<Scalar,int> lhs,
    std::span<const Scalar> solution,
    std::span<const Scalar> rhs,
    Ref<const Mesh> rMesh,
    std::span<const CellData> contiguousCellData,
    Ref<const BVH> rBVH,
    Ref<const Assembler> rAssembler,
    Ref<const utils::ArgParse::Results> rArguments,
    Ref<mp::ThreadPoolBase> rThreads);


} // namespace cie::fem
