#ifndef CIE_FEM_NUMERIC_GAUSS_LEGENDRE_QUADRATURE_HPP
#define CIE_FEM_NUMERIC_GAUSS_LEGENDRE_QUADRATURE_HPP

// --- Internal Includes ---
#include "packages/numeric/inc/QuadratureBase.hpp"


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

    GaussLegendreQuadrature(Size integrationOrder,
                            NT maxAbsoluteNodeError  = 1e-14,
                            Size maxNewtonIterations = 25ul);
}; // class GaussLegendreQuadrature

///@}

} // namespace cie::fem

#endif
