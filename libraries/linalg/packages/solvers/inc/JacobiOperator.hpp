#pragma once

// --- Linalg Includes ---
#include "packages/solvers/inc/LinalgSpace.hpp"
#include "packages/solvers/inc/LinearOperator.hpp"
#include "packages/solvers/inc/DiagonalOperator.hpp"

// --- STL Includes ---
#include <memory>


namespace cie::linalg {


template <LinalgSpaceLike TSpace>
class JacobiOperator
    : public LinearOperator<TSpace> {
public:
    constexpr JacobiOperator() noexcept = default;

    JacobiOperator(
        std::shared_ptr<TSpace> pSpace,
        std::size_t size,
        std::shared_ptr<LinearOperator<TSpace>> pLhs,
        std::shared_ptr<LinearOperator<TSpace>> pInverseDiagonal,
        std::size_t iterations = 1,
        typename TSpace::Value relaxation = 1);

    /// @copydoc LinearOperator::product
    void product(
        typename TSpace::Value inScale,
        typename TSpace::ConstVectorView in,
        typename TSpace::Value outScale,
        typename TSpace::VectorView out) override;

protected:
    std::shared_ptr<LinearOperator<TSpace>> _pLhs, _pInverseDiagonal;

    std::size_t _iterations;

    typename TSpace::Value _relaxation;

    std::shared_ptr<TSpace> _pSpace;

    typename TSpace::Vector _previous, _residual, _memory;
}; // class JacobiOperator


} // namespace cie::linalg
