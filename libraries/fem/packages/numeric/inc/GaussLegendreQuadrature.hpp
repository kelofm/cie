#pragma once

// --- Utility Includes ---
#include "packages/maths/inc/Comparison.hpp"

// --- FEM Includes ---
#include "packages/numeric/inc/QuadratureBase.hpp"

// --- STL Includes ---
#include <limits>


namespace cie::fem {


///@addtogroup fem
///@{

template <concepts::Numeric NT>
class GaussLegendreQuadrature final : public QuadratureBase<NT>
{
public:
    using typename QuadratureBase<NT>::NodeContainer;

    using typename QuadratureBase<NT>::WeightContainer;

public:
    GaussLegendreQuadrature() noexcept = default;

    GaussLegendreQuadrature(
        Size integrationOrder,
        utils::Comparison<NT> comparison = {
            0x100 * std::numeric_limits<NT>::epsilon(),
            0x100 * std::numeric_limits<NT>::epsilon()},
        Size maxNewtonIterations = 50ul);
}; // class GaussLegendreQuadrature

///@}

} // namespace cie::fem
