#pragma once

// help the language server
#include "packages/utilities/inc/ParametricSpace.hpp"

// --- Utility Includes ---
#include "packages/maths/inc/OuterProduct.hpp"

// --- STL Includes ---
#include <array>


namespace cie::fem {


template <unsigned Dimension, concepts::Numeric TValue, ParametricSpaceType ST>
template <concepts::FunctionWithSignature<bool,Ref<const std::span<const std::uint8_t,intPow(2u,Dimension)>>> TFunctor>
constexpr void ParametricSpace<Dimension,TValue,ST>::iterateCorners(TFunctor &&rFunctor) noexcept {
    constexpr unsigned cornerCount = intPow(2u,Dimension);
    std::array<std::uint8_t,cornerCount> state;
    std::fill_n(
        state.data(),
        cornerCount,
        false);
    do {
        if (!rFunctor(state)) break;
    } while (cie::maths::OuterProduct<Dimension>::next(2, state.data()));
}


template <unsigned Dimension, concepts::Numeric TValue, ParametricSpaceType ST>
template <concepts::FunctionWithSignature<bool,Ref<const std::span<const TValue,intPow(2u,Dimension)>>> TFunctor>
constexpr void ParametricSpace<Dimension,TValue,ST>::iterateCorners(TFunctor &&rFunctor) noexcept {
    constexpr unsigned cornerCount = intPow(2u,Dimension);
    ParametricSpace::iterateCorners(
        [&rFunctor] (Ref<const std::span<const std::uint8_t,cornerCount>> state) -> bool {
            std::array<TValue,cornerCount> corner;
            std::transform(
                state.begin(),
                state.end(),
                corner.begin(),
                [](std::uint8_t state) -> TValue {
                    return state ? static_cast<TValue>(1) : static_cast<TValue>(-1);
                });
            return rFunctor(corner);
        });
}


} // namespace cie::fem
