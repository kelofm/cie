#pragma once

// --- Internal Includes ---
#include "poisson2D/mesh.hpp"

// --- FEM Includes ---
#include "packages/graph/inc/Assembler.hpp"


namespace cie::fem {


void integrateStiffness(
    Ref<const Mesh> rMesh,
    std::span<const CellData> contiguousCellData,
    Ref<const Assembler> rAssembler,
    linalg::CSRView<Scalar,int> lhs,
    Ref<const utils::ArgParse::Results> rArguments,
    Ref<mp::ThreadPoolBase> rThreads);


} // namespace cie::fem
