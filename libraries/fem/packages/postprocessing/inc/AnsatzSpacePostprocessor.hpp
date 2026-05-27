#pragma once

// --- FEM Includes ---
#include "packages/numeric/inc/MeshBase.hpp"
#include "packages/graph/inc/Assembler.hpp"

// --- GEO Includes ---
#include "packages/partitioning/inc/AABBoxNode.hpp"

// --- Utility Includes ---
//#include "packages/concurrency/inc/ThreadPoolBase.hpp"

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
        std::span<const ParametricCoordinate<Value>,ParametricDimension> parametricPoint,
        std::span<const Value> fieldValues,
        std::uint8_t fieldComponentCount,
        std::span<const std::size_t> dofIndices,
        std::span<const std::uint8_t> dofOrders,
        std::span<const std::uint8_t> outputOrders,
        std::span<Value> ansatzValueBuffer,
        std::span<Value> ansatzBuffer,
        std::span<Value> out) noexcept;

    template <DiscretizationLike TMesh>
    void interpolate(
        std::span<const PhysicalCoordinate<Value>> physicalPoints,
        Ref<const TMesh> rMesh,
        Ref<const Assembler> rAssembler,
        Ref<const geo::FlatAABBoxTree<Value,PhysicalDimension>> rBVH,
        std::span<const std::span<const Value>> fieldValues,
        std::uint8_t fieldComponentCount,
        std::span<const std::uint8_t> dofOrders,
        std::span<const std::uint8_t> outputOrders,
        std::span<const std::span<Value>> out);

private:
    Value _nanReplacement;
}; // class AnsatzSpacePostprocessor


} // namespace cie::fem

#include "packages/postprocessing/impl/AnsatzSpacePostprocessor_impl.hpp"
