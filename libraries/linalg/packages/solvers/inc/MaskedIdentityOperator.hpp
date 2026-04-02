#pragma once

// --- Linalg Includes ---
#include "packages/solvers/inc/LinearOperator.hpp"
#include "packages/solvers/inc/LinalgSpace.hpp"

// --- STL Includes ---
#include <memory>


namespace cie::linalg {


template <LinalgSpaceLike TSpace>
class MaskedIdentityOperator : public LinearOperator<TSpace> {
public:
    constexpr MaskedIdentityOperator() noexcept = default;

    MaskedIdentityOperator(
        std::shared_ptr<TSpace> pSpace,
        typename TSpace::ConstVectorView mask,
        typename TSpace::Value threshold);

    void product(
        typename TSpace::Value inScale,
        typename TSpace::ConstVectorView in,
        typename TSpace::Value outScale,
        typename TSpace::VectorView out) override;

private:
    typename TSpace::ConstVectorView _mask;

    std::shared_ptr<TSpace> _pSpace;

    typename TSpace::Value _threshold;
}; // class MaskedIdentityOperator


} // namespace cie::linalg
