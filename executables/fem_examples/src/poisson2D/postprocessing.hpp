#pragma once

// --- Internal Includes ---
#include "poisson2D/mesh.hpp"
#include "poisson2D/constraints.hpp"

// --- Utility Includes ---
#include "packages/io/inc/json.hpp"


namespace cie::fem {


void postprocess(
    std::span<const Scalar> meshBase,
    std::span<const Scalar> meshLengths,
    linalg::CSRView<Scalar,int> lhs,
    std::span<const Scalar> solution,
    std::span<const Scalar> rhs,
    Ref<const Mesh> rMesh,
    Ref<const BVH> rBVH,
    std::span<const BoundarySegment> boundarySegments,
    Ref<const Assembler> rAssembler,
    Ref<const cie::io::JSONObject> rConfiguration,
    Ref<mp::ThreadPoolBase> rThreads);


} // namespace cie::fem
