#pragma once

// --- Internal Includes ---
#include "poisson2D/mesh.hpp"

// --- FEM Includes ---
#include "packages/graph/inc/Assembler.hpp"


namespace cie::fem {


void integrateStiffness(
    Ref<const Mesh> rMesh,
    std::span<const Cell> cells,
    Ref<const Assembler> rAssembler,
    linalg::CSRView<Scalar,int> lhs,
    Ref<const cie::io::JSONObject> rConfiguration,
    Ref<mp::ThreadPoolBase> rThreads);


} // namespace cie::fem
