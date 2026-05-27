#pragma once

// --- Internal Includes ---
#include "poisson2D/mesh.hpp"
#include "poisson2D/integration.hpp"

// --- FEM Includes ---
#include "packages/graph/inc/Assembler.hpp"
#include "packages/maths/inc/AffineEmbedding.hpp"
#include "packages/numeric/inc/Quadrature.hpp"
#include "packages/numeric/inc/QuadraturePointFactory.hpp"

// --- GEO Includes ---
#include "packages/partitioning/inc/AABBoxNode.hpp"

// --- Utility Includes ---
#include "packages/io/inc/json.hpp"


namespace cie::fem {


/// @brief Data structure unique to the triangulated, immersed boundary cells.
class BoundaryCell : public CellBase<1,Scalar,maths::AffineEmbedding<Scalar,1u,Dimension>,void,Dimension> {
public:
    using Base = CellBase<1,Scalar,maths::AffineEmbedding<Scalar,1u,Dimension>,void,Dimension>;

    using Base::Base;

    BoundaryCell(
        unsigned id,
        std::size_t ansatzID,
        RightRef<typename Base::SpatialTransform> rEmbedding,
        Ref<const std::array<Scalar,2>> state,
        Ref<const Cell> rCell) noexcept;

    [[nodiscard]] Ref<const Cell> getEmbeddingCell() const noexcept {
        return *_pCell;
    }

    [[nodiscard]] constexpr std::span<const Scalar,2> state() const noexcept {
        return std::span<const Scalar,2>(_state.data(), 2);
    }

private:
    Ptr<const Cell> _pCell;

    std::array<Scalar,2> _state;
}; // class BoundaryCell


class BoundaryMeshData {
public:
    BoundaryMeshData() noexcept = default;

    BoundaryMeshData(
        std::span<const QuadraturePoint<1,Scalar>> quadraturePointSet,
        RightRef<std::vector<BoundaryCell>> rBoundaryCells);

    [[nodiscard]] CachedQuadraturePointFactory<1,Scalar> makeQuadratureRule() const;

    [[nodiscard]] constexpr std::span<const BoundaryCell> cells() const noexcept {
        return _cells;
    }

private:
    std::vector<QuadraturePoint<1,Scalar>> _quadraturePointSet;

    std::vector<BoundaryCell> _cells;
}; // class BoundaryMeshData


/// @brief Mesh type of the immersed, triangulated boundary.
/// @details Cell data consists of
using BoundaryMesh = Graph<
    void,
    void,
    BoundaryMeshData>;


using BVH = geo::FlatAABBoxTree<Scalar,Dimension>;


[[nodiscard]] BVH makeBoundingVolumeHierarchy(
    std::span<Cell> cells,
    std::span<const Scalar> meshBase,
    std::span<const Scalar> meshLengths);


using BoundarySegment = StaticArray<Scalar,2*Dimension+2>;


[[nodiscard]] std::vector<BoundarySegment> makeBoundary(Ref<const cie::io::JSONObject> configuration);


[[nodiscard]] DynamicArray<BoundarySegment>
integrateBoundaryConstraints(
    Ref<Mesh> rMesh,
    std::span<const BoundarySegment> tesselatedBoundary,
    Ref<const Assembler> rAssembler,
    BVH::View bvh,
    std::span<const Cell> cells,
    linalg::CSRView<const Scalar,const int> lhs,
    Ref<DynamicArray<int>> rConstraintRowExtents,
    Ref<DynamicArray<int>> rConstraintColumnIndices,
    Ref<DynamicArray<Scalar>> rConstraintEntries,
    Ref<DynamicArray<Scalar>> rConstraintRHS,
    Ref<const cie::io::JSONObject> rConfiguration);


} // namespace cie::fem
