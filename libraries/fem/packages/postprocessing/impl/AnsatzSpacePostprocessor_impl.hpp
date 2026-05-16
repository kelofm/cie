#pragma once

// --- FEM Includes ---
#include "packages/postprocessing/inc/AnsatzSpacePostprocessor.hpp"

// --- STL Includes ---
#include <limits>


namespace cie::fem {


template <maths::Expression TA, unsigned PD>
AnsatzSpacePostprocessor<TA,PD>::AnsatzSpacePostprocessor() noexcept
    : AnsatzSpacePostprocessor(std::numeric_limits<Value>::quiet_NaN())
{}


template <maths::Expression TA, unsigned PD>
AnsatzSpacePostprocessor<TA,PD>::AnsatzSpacePostprocessor(Value nanReplacement) noexcept
    : _nanReplacement(nanReplacement)
{}


template <maths::Expression TA, unsigned PD>
void AnsatzSpacePostprocessor<TA,PD>::interpolate(
    Ref<const TA> rAnsatzSpace,
    std::span<const Value,ParametricDimension> parametricPoint,
    std::span<const Value> fieldValues,
    std::uint8_t fieldComponentCount,
    std::span<const std::size_t> dofIndices,
    std::span<const std::uint8_t> dofOrders,
    std::span<const std::uint8_t> outputOrders,
    std::span<Value> ansatzValueBuffer,
    std::span<Value> ansatzBuffer,
    std::span<Value> out) noexcept {
        assert(rAnsatzSpace.size() == ansatzValueBuffer.size());
        assert(rAnsatzSpace.bufferSize() <= ansatzBuffer.size());
        assert(dofOrders.size() == rAnsatzSpace.size());
        assert(fieldComponentCount);
        assert(dofIndices.size() % fieldComponentCount == 0);
        assert(dofIndices.size() / fieldComponentCount == rAnsatzSpace.size());
        assert(out.size() == outputOrders.size() * fieldComponentCount);

        rAnsatzSpace.evaluate(
            parametricPoint,
            ansatzValueBuffer,
            ansatzBuffer);

        Ptr<Value> itOut = out.data();
        for (std::uint8_t outputOrder : outputOrders) {
            for (std::uint8_t iFieldComponent=0; iFieldComponent<fieldComponentCount; ++iFieldComponent) {
                for (std::size_t iAnsatz=0ul; iAnsatz<ansatzValueBuffer.size(); ++iAnsatz) {
                    if (dofOrders[iAnsatz] == outputOrder) {
                        (*itOut) += ansatzValueBuffer[iAnsatz] * fieldValues[
                            dofIndices[
                                iFieldComponent * ansatzValueBuffer.size() + iAnsatz]];
                    }
                } // for iAnsatz in range(ansatzValueBuffer.size())
            } // for iFieldComponent in range(fieldComponentCount)
            ++itOut;
        }

        if (_nanReplacement != std::numeric_limits<Value>::quiet_NaN())
            for (Ref<Value> rOut : out)
                if (rOut == std::numeric_limits<Value>::quiet_NaN())
                    rOut = _nanReplacement;
}


template <maths::Expression TA, unsigned PD>
template <CellLike TC, DiscretizationLike TM>
void AnsatzSpacePostprocessor<TA,PD>::interpolate(
    std::span<const Value> parametricPoints,
    std::span<const TC> cells,
    Ref<const TM> rMesh,
    Ref<const geo::FlatAABBoxTree<Value,PD>> rBVH,
    std::span<const Value> fieldValues,
    std::uint8_t fieldComponentCount,
    std::span<Value> out) {
        CIE_BEGIN_EXCEPTION_TRACING

        CIE_END_EXCEPTION_TRACING
}


} // namespace cie::fem
