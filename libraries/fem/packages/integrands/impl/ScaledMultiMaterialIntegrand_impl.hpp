#pragma once

// --- FEM Includes ---
#include "packages/integrands/inc/ScaledMultiMaterialIntegrand.hpp"

// --- Utility Includes ---
#include "packages/io/inc/Serializer.hpp"

// --- STL Includes ---
#include <algorithm>
#include <cassert>


namespace cie::fem {


template <maths::SliceableExpression TI, class MID, int MC>
ScaledMultiMaterialIntegrand<TI,MID,MC>::ScaledMultiMaterialIntegrand(
    RightRef<TI> rIntegrand,
    std::span<const std::pair<MID,Value>,MC> materialMap)
        :   _integrand(std::move(rIntegrand)),
            _materialMap() {
                assert(std::is_sorted(
                    materialMap.begin(),
                    materialMap.end(),
                    [] (const auto& rLeft, const auto& rRight) {
                        return rLeft.first < rRight.first;
                    }));
                std::copy_n(
                    materialMap.data(),
                    MC,
                    _materialMap.data());
}


/// @name Common Expression Interface
/// @{


template <maths::SliceableExpression TI, class MID, int MC>
void ScaledMultiMaterialIntegrand<TI,MID,MC>::evaluate(
    ConstSpan in,
    Span out,
    BufferSpan buffer) const {
        const auto bounds = this->contextBounds();
        Context context;
        for (unsigned iDimension=0u; iDimension<ParametricDimension; ++iDimension) {
            context.data()[2 * iDimension] = 0;
            context.data()[2 * iDimension + 1] = bounds[iDimension];
        }
        InContext inContext;
        inContext.fill(context);
        return this->evaluate(in, out, buffer, inContext);
}


template <maths::SliceableExpression TI, class MID, int MC>
void ScaledMultiMaterialIntegrand<TI,MID,MC>::evaluate(
    ConstSpan in,
    Span out,
    BufferSpan buffer,
    Ref<const Context> rContext) const {
        InContext inContext;
        inContext.fill(rContext);
        return this->evaluate(in, out, buffer, inContext);
}


template <maths::SliceableExpression TI, class MID, int MC>
void ScaledMultiMaterialIntegrand<TI,MID,MC>::evaluate(
    ConstSpan in,
    Span out,
    BufferSpan buffer,
    Ref<const InContext> rInContext) const {
        // Assume the material ID is the last entry in the input.
        Value scale = 0;
        const MID materialID = static_cast<MID>(in.back());
        const auto itMaterialSentinel = std::upper_bound(
            _materialMap.begin(),
            _materialMap.end(),
            materialID,
            [] (MID materialID, const auto& rPair) {return materialID < rPair.first;});
        if (itMaterialSentinel == _materialMap.begin()) {
            scale = _materialMap.front().second;
        } else if (itMaterialSentinel == _materialMap.end()) {
            scale = _materialMap.back().second;
        } else if constexpr (1 < MC) {
            const auto itMaterial = itMaterialSentinel - 1;
            scale = itMaterial->second + (itMaterialSentinel->second - itMaterial->second) / (itMaterialSentinel->first - itMaterial->first) * (materialID - itMaterial->first);
        } else static_assert(MC == 1, "invalid material map");
        _integrand.evaluate(
            {in.data(), in.size() - 1},
            out,
            buffer,
            rInContext);
        for (Ref<Value> v : out) v *= scale;
}


/// @}
/// @name Runtime Expression Interface
/// @{


template <maths::SliceableExpression TI, class MID, int MC>
unsigned ScaledMultiMaterialIntegrand<TI,MID,MC>::size() const noexcept
requires (!maths::StaticExpression<TI>) {
    return _integrand.size();
}


template <maths::SliceableExpression TI, class MID, int MC>
unsigned ScaledMultiMaterialIntegrand<TI,MID,MC>::bufferSize() const noexcept
requires (!maths::StaticExpression<TI>) {
    return _integrand.bufferSize();
}


template <maths::SliceableExpression TI, class MID, int MC>
unsigned ScaledMultiMaterialIntegrand<TI,MID,MC>::size(Ref<const Context> rContext) const noexcept
requires (!maths::StaticExpression<TI>) {
    InContext inContext;
    inContext.fill(rContext);
    return this->size(inContext);
}


template <maths::SliceableExpression TI, class MID, int MC>
unsigned ScaledMultiMaterialIntegrand<TI,MID,MC>::size(Ref<const InContext> rInContext) const noexcept
requires (!maths::StaticExpression<TI>) {
    return _integrand.size(rInContext);
}


template <maths::SliceableExpression TI, class MID, int MC>
unsigned ScaledMultiMaterialIntegrand<TI,MID,MC>::bufferSize(Ref<const Context> rContext) const noexcept
requires (!maths::StaticExpression<TI>) {
    InContext inContext;
    inContext.fill(rContext);
    return this->bufferSize(inContext);
}


template <maths::SliceableExpression TI, class MID, int MC>
unsigned ScaledMultiMaterialIntegrand<TI,MID,MC>::bufferSize(Ref<const InContext> rInContext) const noexcept
requires (!maths::StaticExpression<TI>) {
    return _integrand.bufferSize(rInContext);
}


template <maths::SliceableExpression TI, class MID, int MC>
std::array<typename TI::Context::Value,ScaledMultiMaterialIntegrand<TI,MID,MC>::ParametricDimension>
ScaledMultiMaterialIntegrand<TI,MID,MC>::contextBounds() const noexcept
requires (!maths::StaticExpression<TI>) {
    return _integrand.contextBounds();
}


template <maths::SliceableExpression TI, class MID, int MC>
typename ScaledMultiMaterialIntegrand<TI,MID,MC>::OutContext
ScaledMultiMaterialIntegrand<TI,MID,MC>::outContext(Ref<const InContext> rInContext) const noexcept
requires (!maths::StaticExpression<TI>) {
    return _integrand.outContext(rInContext);
}


/// @}
/// @name Static Expression Interface
/// @{


template <maths::SliceableExpression TI, class MID, int MC>
constexpr unsigned ScaledMultiMaterialIntegrand<TI,MID,MC>::size() noexcept
requires maths::StaticExpression<TI> {
    return TI::size();
}


template <maths::SliceableExpression TI, class MID, int MC>
constexpr unsigned ScaledMultiMaterialIntegrand<TI,MID,MC>::bufferSize() noexcept
requires maths::StaticExpression<TI> {
    return TI::bufferSize();
}


template <maths::SliceableExpression TI, class MID, int MC>
constexpr unsigned ScaledMultiMaterialIntegrand<TI,MID,MC>::size(Ref<const Context> rContext) noexcept
requires maths::StaticExpression<TI> {
    InContext inContext;
    inContext.fill(rContext);
    return size(inContext);
}


template <maths::SliceableExpression TI, class MID, int MC>
constexpr unsigned ScaledMultiMaterialIntegrand<TI,MID,MC>::size(Ref<const InContext> rInContext) noexcept
requires maths::StaticExpression<TI> {
    return TI::size(rInContext);
}


template <maths::SliceableExpression TI, class MID, int MC>
constexpr unsigned ScaledMultiMaterialIntegrand<TI,MID,MC>::bufferSize(Ref<const Context> rContext) noexcept
requires maths::StaticExpression<TI> {
    InContext inContext;
    inContext.fill(rContext);
    return bufferSize(inContext);
}


template <maths::SliceableExpression TI, class MID, int MC>
constexpr unsigned ScaledMultiMaterialIntegrand<TI,MID,MC>::bufferSize(Ref<const InContext> rInContext) noexcept
requires maths::StaticExpression<TI> {
    return TI::bufferSize(rInContext);
}


template <maths::SliceableExpression TI, class MID, int MC>
constexpr std::array<typename TI::Context::Value,ScaledMultiMaterialIntegrand<TI,MID,MC>::ParametricDimension>
ScaledMultiMaterialIntegrand<TI,MID,MC>::contextBounds() noexcept
requires (maths::StaticExpression<TI>) {
    return TI::contextBounds();
}


template <maths::SliceableExpression TI, class MID, int MC>
constexpr typename ScaledMultiMaterialIntegrand<TI,MID,MC>::OutContext
ScaledMultiMaterialIntegrand<TI,MID,MC>::outContext(Ref<const InContext> rInContext) noexcept
requires maths::StaticExpression<TI> {
    return TI::outContext(rInContext);
}


/// @}
/// @name Serialization
/// @{


template <maths::SliceableExpression TI, class MID, int MC>
void ScaledMultiMaterialIntegrand<TI,MID,MC>::serialize(
    Ref<cie::io::Traits::SerializerStream> rStream,
    tags::Binary) const {
        using BS = cie::io::BinarySerializer;
        BS::serialize(
            rStream,
            _integrand);
        BS::serialize(
            rStream,
            _materialMap.data(),
            _materialMap.size());
}


template <maths::SliceableExpression TI, class MID, int MC>
void ScaledMultiMaterialIntegrand<TI,MID,MC>::deserialize(
    Ref<cie::io::Traits::DeserializerStream> rStream,
    Ref<ScaledMultiMaterialIntegrand> rInstance,
    tags::Binary) {
        using BS = cie::io::BinarySerializer;
        BS::deserialize(
            rStream,
            rInstance._integrand);
        BS::deserialize(
            rStream,
            rInstance._materialMap.data(),
            rInstance._materialMap.size());
}


/// @}


} // namespace cie::fem
