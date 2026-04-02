#pragma once

// --- Linalg Includes ---
#include "packages/solvers/inc/LinalgSpace.hpp"


namespace cie::linalg {


template <LinalgSpaceLike TSpace>
class LinearOperator {
public:
    virtual ~LinearOperator() = default;

    /// @brief Compute a scaled matrix-vector product and add it to the provided output vector.
    /// @details Computes @f[
    ///             r = \alpha r + \beta A b
    ///          @f]
    ///          where
    ///          - @p A is the linear operator represented by this class instance,
    ///          - @p b is the input vector to compute the product of,
    ///          - @p r is the output vector to add the scaled product to.
    ///          - @p @f$\alpha@f$ is the scaling factor for the initial value of the output vector @f$r@f$,
    ///          - @p @f$\beta@f$ is the scaling factor for the matrix-vector product @f$A b@f$.
    /// @param inScale Factor to scale the initial value of the output vector with (@f$\alpha@f$).
    /// @param in Input vector to compute the product of (@f$b@f$).
    /// @param outScale Factor to scale the matrix-vector product with (@f$\beta@f$).
    /// @param out Output array to add the scaled product to (@f$r@f$).
    virtual void product(
        typename TSpace::Value inScale,
        typename TSpace::ConstVectorView in,
        typename TSpace::Value outScale,
        typename TSpace::VectorView out) = 0;
}; // class LinearOperator


} // namespace cie::linalg
