#pragma once

// --- Internal Includes ---
#include "poisson2D/mesh.hpp"

// --- FEM Includes ---
#include "packages/graph/inc/Assembler.hpp"


namespace cie::fem {


struct CSRWrapper {
    const int rowCount, columnCount;
    std::span<const int> rowExtents, columnIndices;
    std::span<Scalar> entries;
}; // struct CSRWrapper


void integrateStiffness(
    Ref<const Mesh> rMesh,
    std::span<const CellData> contiguousCellData,
    Ref<const Assembler> rAssembler,
    CSRWrapper lhs,
    Ref<const utils::ArgParse::Results> rArguments,
    Ref<mp::ThreadPoolBase> rThreads);


} // namespace cie::fem
