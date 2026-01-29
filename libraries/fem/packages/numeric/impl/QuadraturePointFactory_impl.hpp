#pragma once

// help the language server
#include "packages/numeric/inc/QuadraturePointFactory.hpp"

// --- Utility Includes ---
#include "packages/maths/inc/OuterProduct.hpp"

// --- STL Includes ---
#include <algorithm>


namespace cie::fem {


template <unsigned Dim, concepts::Numeric TValue, class TQD>
constexpr inline OuterProductQuadraturePointFactory<Dim,TValue,TQD>::OuterProductQuadraturePointFactory() noexcept
    : OuterProductQuadraturePointFactory(std::span<const QuadraturePoint<1u,TValue,TQD>>())
{}


template <unsigned Dim, concepts::Numeric TValue, class TQD>
constexpr inline OuterProductQuadraturePointFactory<Dim,TValue,TQD>::OuterProductQuadraturePointFactory(std::span<const QuadraturePoint<1u,TValue,TQD>> base) noexcept
    : _done(base.empty()),
      _base(base),
      _state()
{
    std::fill_n(
        _state.data(),
        Dimension,
        0u);
}


template <unsigned Dim, concepts::Numeric TValue, class TQD>
unsigned OuterProductQuadraturePointFactory<Dim,TValue,TQD>::operator()(std::span<Value> out) noexcept {
    auto itOut = out.begin();

    // Early exit if all quadrature points were already generated.
    if (_done || itOut == out.end()) return 0u;

    // Generate quadrature points until they're exhausted
    // or there's no more output slots.
    do {
        Ref<Value> rCurrent = *itOut++;
        rCurrent.weight() = static_cast<TValue>(1);
        for (unsigned iState=0u; iState<Dimension; ++iState) {
            rCurrent.weight() *= _base[_state[iState]].weight();
            rCurrent.position()[iState] = _base[_state[iState]].position().front();
        }
    } while (cie::maths::OuterProduct<Dimension>::next(_base.size(), _state.data()) && itOut != out.end());

    // Check whether all quadrature points have been generated.
    _done = std::all_of(
        _state.begin(),
        _state.end(),
        [](unsigned iNode){return iNode == 0u;});

    return std::distance(out.begin(), itOut);
}


template <unsigned Dim, concepts::Numeric TValue, class TQD>
CachedQuadraturePointFactory<Dim,TValue,TQD>::CachedQuadraturePointFactory() noexcept
    : CachedQuadraturePointFactory(std::span<const Value>())
{}


template <unsigned Dim, concepts::Numeric TValue, class TQD>
CachedQuadraturePointFactory<Dim,TValue,TQD>::CachedQuadraturePointFactory(std::span<const Value> cache) noexcept
    : _cache(cache)
{}


template <unsigned Dim, concepts::Numeric TValue, class TQD>
unsigned CachedQuadraturePointFactory<Dim,TValue,TQD>::operator()(std::span<Value> out) noexcept {
    const std::size_t copyCount = std::min(_cache.size(), out.size());
    std::copy_n(
        _cache.data(),
        copyCount,
        out.data());
    _cache = std::span<const Value>(
        _cache.data() + copyCount,
        _cache.size() - copyCount);
    return copyCount;
}


} // namespace cie::fem
