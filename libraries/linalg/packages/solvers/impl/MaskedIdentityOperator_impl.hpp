#pragma once

// --- Linalg Includes ---
#include "packages/solvers/inc/MaskedIdentityOperator.hpp"


namespace cie::linalg {


template <LinalgSpaceLike TSpace>
MaskedIdentityOperator<TSpace>::MaskedIdentityOperator(
    std::shared_ptr<TSpace> pSpace,
    typename TSpace::ConstVectorView mask,
    typename TSpace::Value threshold)
        :   _mask(mask),
            _pSpace(pSpace),
            _threshold(threshold)
{}


} // namespace cie::linalg
