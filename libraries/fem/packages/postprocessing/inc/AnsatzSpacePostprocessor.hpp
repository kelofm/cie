#pragma once

// --- FEM Includes ---
#include "packages/maths/inc/AnsatzSpace.hpp"
#include "packages/maths/inc/Polynomial.hpp"
#include "packages/numeric/inc/CellBase.hpp"
#include "packages/numeric/inc/MeshBase.hpp"

// --- GEO Includes ---
#include "packages/partitioning/inc/AABBoxNode.hpp"

// --- Utility Includes ---
#include "packages/concurrency/inc/ThreadPoolBase.hpp"

// --- STL Includes ---
#include <span>


namespace cie::fem {


template <
    maths::Expression TAnsatz,
    unsigned PhysicalDimension = TAnsatz::Dimension>
class AnsatzSpacePostprocessor {
public:
    using Ansatz = TAnsatz;

    using Value = typename TAnsatz::Value;

    static constexpr inline unsigned ParametricDimension = Ansatz::Dimension;

    AnsatzSpacePostprocessor() noexcept;

    AnsatzSpacePostprocessor(Value nanReplacement) noexcept;

    void interpolate(
        Ref<const TAnsatz> rAnsatzSpace,
        std::span<const Value,ParametricDimension> parametricPoint,
        std::span<const Value> fieldValues,
        std::uint8_t fieldComponentCount,
        std::span<const std::size_t> dofIndices,
        std::span<const std::uint8_t> dofOrders,
        std::span<const std::uint8_t> outputOrders,
        std::span<Value> ansatzValueBuffer,
        std::span<Value> ansatzBuffer,
        std::span<Value> out) noexcept;

    template <CellLike TCell, DiscretizationLike TMesh>
    void interpolate(
        std::span<const Value> parametricPoints,
        std::span<const TCell> cells,
        Ref<const TMesh> rMesh,
        Ref<const geo::FlatAABBoxTree<Value,PhysicalDimension>> rBVH,
        std::span<const Value> fieldValues,
        std::uint8_t fieldComponentCount,
        std::span<Value> out);

private:
    Value _nanReplacement;
}; // class AnsatzSpacePostprocessor


} // namespace cie::fem
