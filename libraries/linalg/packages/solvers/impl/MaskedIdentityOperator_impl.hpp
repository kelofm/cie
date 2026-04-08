#pragma once

// --- Linalg Includes ---
#include "packages/solvers/inc/MaskedIdentityOperator.hpp"


namespace cie::linalg {


template <LinalgSpaceLike TSpace, LinalgSpaceLike TMaskSpace>
MaskedIdentityOperator<TSpace,TMaskSpace>::MaskedIdentityOperator(
    std::shared_ptr<TSpace> pSpace,
    std::shared_ptr<TMaskSpace> pMaskSpace,
    typename TMaskSpace::ConstVectorView mask,
    typename TMaskSpace::Value threshold)
        :   _mask(mask),
            _pSpace(pSpace),
            _pMaskSpace(pMaskSpace),
            _threshold(threshold)
{}


} // namespace cie::linalg
