#pragma once

// --- FEM Includes ---
#include "packages/integrands/inc/ScaledMultiMaterialIntegrand.hpp"

// --- STL Includes ---
#include <algorithm>
#include <cassert>


namespace cie::fem {


template <maths::Expression TI, class MID, int MC>
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


template <maths::Expression TI, class MID, int MC>
void ScaledMultiMaterialIntegrand<TI,MID,MC>::evaluate(
    ConstSpan in,
    Span out,
    BufferSpan buffer) const {
        // Assume the material ID is the last entry in the input.
        assert(in.size() == _integrand.size() + 1);
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
            buffer);
        for (Ref<Value> v : out) v *= scale;
}


template <maths::Expression TI, class MID, int MC>
unsigned ScaledMultiMaterialIntegrand<TI,MID,MC>::size() const noexcept
requires (!maths::StaticExpression<TI>) {
    return _integrand.size();
}


template <maths::Expression TI, class MID, int MC>
constexpr unsigned ScaledMultiMaterialIntegrand<TI,MID,MC>::size() noexcept
requires maths::StaticExpression<TI> {
    return TI::size();
}


template <maths::Expression TI, class MID, int MC>
unsigned ScaledMultiMaterialIntegrand<TI,MID,MC>::bufferSize() const noexcept
requires (!maths::StaticExpression<TI>) {
    return _integrand.bufferSize();
}


template <maths::Expression TI, class MID, int MC>
constexpr unsigned ScaledMultiMaterialIntegrand<TI,MID,MC>::bufferSize() noexcept
requires maths::StaticExpression<TI> {
    return TI::bufferSize();
}


} // namespace cie::fem
