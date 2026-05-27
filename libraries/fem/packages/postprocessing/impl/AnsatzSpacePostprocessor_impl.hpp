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
    std::span<const ParametricCoordinate<Value>,ParametricDimension> parametricPoint,
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
        std::fill(
            out.begin(),
            out.end(),
            Value(0));

        rAnsatzSpace.evaluate(
            Kernel<PD,Value>::decay(parametricPoint),
            ansatzValueBuffer,
            ansatzBuffer);

        for (std::uint8_t iFieldComponent=0; iFieldComponent<fieldComponentCount; ++iFieldComponent) {
            for (std::size_t iOutputOrder=0ul; iOutputOrder<outputOrders.size(); ++iOutputOrder) {
                const auto outputOrder = outputOrders[iOutputOrder];
                for (std::size_t iAnsatz=0ul; iAnsatz<ansatzValueBuffer.size(); ++iAnsatz) {
                    if (dofOrders[iAnsatz] == outputOrder) {
                        out[iFieldComponent * outputOrders.size() + iOutputOrder] += ansatzValueBuffer[iAnsatz] * fieldValues[
                            dofIndices[iFieldComponent * ansatzValueBuffer.size() + iAnsatz]];
                    }
                } // for iAnsatz in range(ansatzValueBuffer.size())
            } // for outputOrder
        } // for iFieldComponent in range(fieldComponentCount)

        if (_nanReplacement != std::numeric_limits<Value>::quiet_NaN())
            for (Ref<Value> rOut : out)
                if (rOut == std::numeric_limits<Value>::quiet_NaN())
                    rOut = _nanReplacement;
}


template <maths::Expression TA, unsigned PhysicalDimension>
template <DiscretizationLike TM>
void AnsatzSpacePostprocessor<TA,PhysicalDimension>::interpolate(
    std::span<const PhysicalCoordinate<Value>> physicalPoints,
    Ref<const TM> rMesh,
    Ref<const Assembler> rAssembler,
    Ref<const geo::FlatAABBoxTree<Value,PhysicalDimension>> rBVH,
    std::span<const std::span<const Value>> fieldValues,
    std::uint8_t fieldComponentCount,
    std::span<const std::uint8_t> dofOrders,
    std::span<const std::uint8_t> outputOrders,
    std::span<const std::span<Value>> out) {
        CIE_CHECK(
            physicalPoints.size() % PhysicalDimension == 0,
            std::format(
                "expecting {}-D points, but got {} components",
                PhysicalDimension,
                physicalPoints.size()));

        CIE_CHECK(
            fieldValues.size() == out.size(),
            std::format(
                "got {} fields, {} output arrays",
                fieldValues.size(),
                out.size()))

        const auto& rCells = rMesh.data().cells();
        const std::size_t pointCount = physicalPoints.size() / PhysicalDimension;
        const auto& bvh = rBVH.makeView();

        //using TCell = typename std::remove_cvref_t<decltype(rCells)>::value_type;
        //constexpr unsigned physicalDimension = TCell::PhysicalDimension;
        for (std::size_t iOutSpan=0ul; iOutSpan<out.size(); ++iOutSpan) {
            const auto span = out[iOutSpan];
            CIE_CHECK(
                span.size() == pointCount * fieldComponentCount * outputOrders.size(),
                std::format(
                    "expecting output array {} of size {} for {} sample points, but got an array of size {}",
                    iOutSpan,
                    pointCount * fieldComponentCount * outputOrders.size(),
                    pointCount,
                    span.size()));
        }
        using PhysicalKernel = Kernel<PhysicalDimension,Value>;

        CIE_BEGIN_EXCEPTION_TRACING
            const auto kernel = [&] (
                std::size_t iPoint,
                Ref<std::vector<Value>> rOutputBuffer,
                Ref<std::vector<std::size_t>> rDoFIndexBuffer,
                Ref<std::vector<Value>> rAnsatzValueBuffer,
                Ref<std::vector<Value>> rAnsatzBuffer) -> void {
                    const std::span<const PhysicalCoordinate<Value>,PhysicalDimension> physicalCoordinates(
                        physicalPoints.data() + iPoint * PhysicalDimension,
                        PhysicalDimension);
                    const std::size_t iMaybeCell = bvh.find(
                        PhysicalKernel::decay(physicalCoordinates),
                        rCells);

                    if (iMaybeCell != rCells.size()) {
                        // Transform the point from physical space to the
                        // found cell's parametric space.
                        const auto& rCell = rCells[iMaybeCell];
                        std::array<ParametricCoordinate<Value>,ParametricDimension> parametricCoordinates;
                        rCell.transform(physicalCoordinates, parametricCoordinates, {});

                        // Resize required buffers.
                        const auto& rDoFs = rAssembler[rCell.id()];
                        Ref<const TA> rAnsatz = rMesh.data().ansatz(rCell.ansatzID());
                        rDoFIndexBuffer.resize(rDoFs.size());
                        rAnsatzValueBuffer.resize(rAnsatz.size());
                        rAnsatzBuffer.resize(rAnsatz.bufferSize());

                        // Collect the cell's DoFs.
                        std::copy(
                            rDoFs.begin(),
                            rDoFs.end(),
                            rDoFIndexBuffer.begin());

                        // Perform the interpolation.
                        std::fill(
                            rOutputBuffer.begin(),
                            rOutputBuffer.end(),
                            static_cast<Value>(0));

                        for (std::size_t iField=0ul; iField<out.size(); ++iField) {
                            const auto fieldSpan = fieldValues[iField];
                            const auto outSpan = out[iField];

                            this->interpolate(
                                rAnsatz,
                                parametricCoordinates,
                                fieldSpan,
                                fieldComponentCount,
                                rDoFIndexBuffer,
                                dofOrders,
                                outputOrders,
                                rAnsatzValueBuffer,
                                rAnsatzBuffer,
                                rOutputBuffer);

                            // Map interpolated values to the output array.
                            // rOutputBuffer[iFieldComponent][iOrder]
                            // outSpan[iFieldComponent][iPoint][iOrder]
                            for (std::size_t iFieldComponent=0ul; iFieldComponent<fieldComponentCount; ++iFieldComponent) {
                                std::span<Value> inRange(
                                    rOutputBuffer.data() + iFieldComponent * outputOrders.size(),
                                    outputOrders.size());
                                std::span<Value> outRange(
                                    outSpan.data() + iFieldComponent * pointCount * outputOrders.size() + iPoint * outputOrders.size(),
                                    outputOrders.size());
                                std::copy(
                                    inRange.begin(),
                                    inRange.end(),
                                    outRange.begin());
                            } // for iFieldComponent in fieldComponentCount
                        } // for iField in out.size()
                    } /*if iMaybeCell != rCells.size()*/ else {
                        for (std::size_t iField=0ul; iField<out.size(); ++iField)
                            for (std::size_t iFieldComponent=0ul; iFieldComponent<fieldComponentCount; ++iFieldComponent) {
                                std::span<Value> outRange(
                                    out[iField].data() + iFieldComponent * pointCount * outputOrders.size() + iPoint * outputOrders.size(),
                                    outputOrders.size());
                                std::fill(
                                    outRange.begin(),
                                    outRange.end(),
                                    _nanReplacement);
                            } // for iFieldComponent
                    }
            }; // kernel

            std::vector<Value>
                outputBuffer(fieldComponentCount * outputOrders.size()),
                ansatzValueBuffer,
                ansatzBuffer;
            std::vector<std::size_t> dofIndexBuffer;
            for (std::size_t iPoint=0ul; iPoint<pointCount; ++iPoint) {
                kernel(
                    iPoint,
                    outputBuffer,
                    dofIndexBuffer,
                    ansatzValueBuffer,
                    ansatzBuffer);
            } // for iPoint
        CIE_END_EXCEPTION_TRACING
}


} // namespace cie::fem
