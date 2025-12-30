#pragma once

// --- FEM Includes ---
#include "packages/utilities/inc/kernel.hpp"

// --- Utility Includes ---
#include "packages/compile_time/packages/concepts/inc/basic_concepts.hpp"
#include "packages/stl_extension/inc/StaticArray.hpp"

// --- STL Includes ---
#include <algorithm>


namespace cie::fem {


template <unsigned Dim, concepts::Numeric TValue>
class QuadraturePoint {
public:
    using Value = typename Kernel<Dim,TValue>::LocalCoordinate;

    constexpr static inline unsigned Dimension = Dim;

    constexpr QuadraturePoint() noexcept {
        std::fill_n(
            _data.data(),
            Dim + 1,
            static_cast<TValue>(0));
    }

    constexpr QuadraturePoint(Ref<const std::span<const TValue,Dim>> rPosition,
                              TValue weight) noexcept {
        std::copy_n(
            rPosition.data(),
            Dim,
            _data.data());
        _data.back() = weight;
    }

    constexpr std::span<const TValue,Dim> position() const noexcept {
        return std::span<const TValue,Dim>(_data.data(), Dim);
    }

    constexpr Ref<const TValue> weight() const noexcept {
        return _data.back();
    }

    constexpr std::span<TValue,Dim> position() noexcept {
        return std::span<TValue,Dim>(_data.data(), Dim);
    }

    constexpr Ref<TValue> weight() noexcept {
        return _data.back();
    }

private:
    StaticArray<TValue,Dim+1> _data;
};


} // namespace cie::fem
