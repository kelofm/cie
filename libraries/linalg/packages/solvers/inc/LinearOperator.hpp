#pragma once

// --- Linalg Includes ---
#include "packages/solvers/inc/LinalgSpace.hpp"


namespace cie::linalg {


template <LinalgSpaceLike TSpace>
class LinearOperator {
public:
    /// @brief Compute a scaled matrix-vector product and add it to the provided output vector.
    /// @details Computes @f[
    ///             r += \alpha * A @ b
    ///          @f]
    ///          where
    ///          - @p A is the linear operator represented by this class instance,
    ///          - @p b is the input vector to compute the product of,
    ///          - @p @f$\alpha@f$ is the constant scaling factor,
    ///          - @p r is the output vector to add the scaled product to.
    /// @param in Input vector to compute the product of (@f$b@f$).
    /// @param scale Factor to scale the product with (@f$\alpha@f$).
    /// @param out Output array to add the scaled product to (@f$r@f$).
    virtual void product(
        typename TSpace::ConstVectorView in,
        typename TSpace::Value scale,
        typename TSpace::VectorView out) const = 0;
}; // class LinearOperator


} // namespace cie::linalg
