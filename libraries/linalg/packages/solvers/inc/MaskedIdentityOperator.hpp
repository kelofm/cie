#pragma once

// --- Linalg Includes ---
#include "packages/solvers/inc/LinearOperator.hpp"
#include "packages/solvers/inc/LinalgSpace.hpp"

// --- STL Includes ---
#include <memory>


namespace cie::linalg {


template <LinalgSpaceLike TSpace, LinalgSpaceLike TMaskSpace>
class MaskedIdentityOperator : public LinearOperator<TSpace> {
public:
    constexpr MaskedIdentityOperator() noexcept = default;

    MaskedIdentityOperator(
        std::shared_ptr<TSpace> pSpace,
        std::shared_ptr<TMaskSpace> pMaskSpace,
        typename TMaskSpace::ConstVectorView mask,
        typename TMaskSpace::Value threshold);

    void product(
        typename TSpace::Value inScale,
        typename TSpace::ConstVectorView in,
        typename TSpace::Value outScale,
        typename TSpace::VectorView out) override;

private:
    typename TMaskSpace::ConstVectorView _mask;

    std::shared_ptr<TSpace> _pSpace;

    std::shared_ptr<TMaskSpace> _pMaskSpace;

    typename TMaskSpace::Value _threshold;
}; // class MaskedIdentityOperator


} // namespace cie::linalg
