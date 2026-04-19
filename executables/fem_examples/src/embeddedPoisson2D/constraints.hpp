#pragma once

// --- Internal Includes ---
#include "embeddedPoisson2D/mesh.hpp"
#include "embeddedPoisson2D/integration.hpp"

// --- FEM Includes ---
#include "packages/graph/inc/Assembler.hpp"
#include "packages/maths/inc/AffineEmbedding.hpp"
#include "packages/numeric/inc/Quadrature.hpp"
#include "packages/numeric/inc/QuadraturePointFactory.hpp"

// --- GEO Includes ---
#include "packages/partitioning/inc/AABBoxNode.hpp"


namespace cie::fem {


/// @brief Data structure unique to the triangulated, immersed boundary cells.
class BoundaryCellData : public CellBase<1,Scalar,maths::AffineEmbedding<Scalar,1u,Dimension>,void,Dimension> {
public:
    using Base = CellBase<1,Scalar,maths::AffineEmbedding<Scalar,1u,Dimension>,void,Dimension>;

    using Base::Base;

    BoundaryCellData(
        unsigned id,
        std::size_t ansatzID,
        RightRef<typename Base::SpatialTransform> rEmbedding,
        Ref<const std::array<Scalar,2>> state,
        Ref<const CellData> rCell) noexcept;

    Ref<const CellData> getEmbeddingCell() const noexcept {
        return *_pCell;
    }

    constexpr std::span<const Scalar,2> state() const noexcept {
        return std::span<const Scalar,2>(_state.data(), 2);
    }

private:
    Ptr<const CellData> _pCell;

    std::array<Scalar,2> _state;
}; // class BoundaryCellData


class BoundaryMeshData {
public:
    BoundaryMeshData() noexcept = default;

    BoundaryMeshData(std::span<const QuadraturePoint<1,Scalar>> quadraturePointSet);

    CachedQuadraturePointFactory<1,Scalar> makeQuadratureRule() const;

private:
    std::vector<QuadraturePoint<1,Scalar>> _quadraturePointSet;
}; // class BoundaryMeshData


/// @brief Mesh type of the immersed, triangulated boundary.
/// @details Cell data consists of
using BoundaryMesh = Graph<
    BoundaryCellData,
    void,
    BoundaryMeshData>;


using BVH = geo::FlatAABBoxTree<Scalar,Dimension>;


BVH makeBoundingVolumeHierarchy(
    std::span<CellData> cells,
    std::span<const Scalar> meshBase,
    std::span<const Scalar> meshLengths);


using BoundarySegment = StaticArray<Scalar,2*Dimension+2>;


std::vector<BoundarySegment> makeBoundary(Ref<const utils::ArgParse::Results> rArguments);


BoundaryMesh generateBoundaryMesh(
    std::span<const BoundarySegment> tesselatedBoundary,
    BVH::View bvh,
    std::span<const CellData> contiguousCellData,
    Ref<const utils::ArgParse::Results> rArguments,
    Ref<DynamicArray<BoundarySegment>> rBoundarySegments);


[[nodiscard]] DynamicArray<BoundarySegment>
imposeBoundaryConditions(
    Ref<Mesh> rMesh,
    std::span<const BoundarySegment> tesselatedBoundary,
    Ref<const Assembler> rAssembler,
    BVH::View bvh,
    std::span<const CellData> contiguousCellData,
    linalg::CSRView<Scalar,int> lhs,
    std::span<Scalar> rhs,
    Ref<const utils::ArgParse::Results> rArguments);


} // namespace cie::fem
