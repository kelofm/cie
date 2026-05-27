#pragma once

// --- Linalg Includes ---
#include "packages/solvers/inc/LinalgSpace.hpp"
#include "packages/solvers/inc/LinearOperator.hpp"

// --- STL Includes ---
#include <memory>


namespace cie::linalg {


template <LinalgSpaceLike TSpace>
class SumOperator : public LinearOperator<TSpace> {
public:
    constexpr SumOperator() noexcept = default;

    SumOperator(
        std::shared_ptr<TSpace> pSpace,
        std::shared_ptr<LinearOperator<TSpace>> pLhs,
        std::shared_ptr<LinearOperator<TSpace>> pRhs,
        std::size_t rowCount);

    void product(
        typename TSpace::Value inScale,
        typename TSpace::ConstVectorView in,
        typename TSpace::Value outScale,
        typename TSpace::VectorView out) override;

private:
    std::shared_ptr<TSpace> _pSpace;

    std::shared_ptr<LinearOperator<TSpace>> _pLhs, _pRhs;

    typename TSpace::Vector _buffer;
}; // class SumOperator


} // namespace linalg

#include "packages/solvers/impl/SumOperator_impl.hpp"
