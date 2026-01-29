#pragma once

// --- Utility Includes ---
#include "packages/maths/inc/power.hpp"
#include "packages/compile_time/packages/concepts/inc/functional.hpp"

// --- STL Includes ---
#include <span>
#include <cstdint>


namespace cie::fem {


enum class ParametricSpaceType {
    Cartesian
}; // enum class ParametricSpaceType


template <
    unsigned Dimension,
    concepts::Numeric TValue,
    ParametricSpaceType SpaceType = ParametricSpaceType::Cartesian>
class ParametricSpace {
public:
    template <concepts::FunctionWithSignature<bool,Ref<const std::span<const std::uint8_t,intPow(2u,Dimension)>>> TFunctor>
    static constexpr void iterateCorners(TFunctor&& rFunctor) noexcept;

    template <concepts::FunctionWithSignature<bool,Ref<const std::span<const TValue,intPow(2u,Dimension)>>> TFunctor>
    static constexpr void iterateCorners(TFunctor&& rFunctor) noexcept;
};


} // namespace cie::fem

#include "packages/utilities/impl/ParametricSpace_impl.hpp"
