#pragma once

// --- Linalg Includes ---
#include "packages/solvers/inc/SumOperator.hpp"


namespace cie::linalg {


template <LinalgSpaceLike TS>
SumOperator<TS>::SumOperator(
    std::shared_ptr<TS> pSpace,
    std::shared_ptr<LinearOperator<TS>> pLhs,
    std::shared_ptr<LinearOperator<TS>> pRhs,
    std::size_t rowCount)
        :   _pSpace(pSpace),
            _pLhs(pLhs),
            _pRhs(pRhs),
            _buffer(pSpace->makeVector(rowCount)) {
                pSpace->fill(_buffer, 0);
}


template <LinalgSpaceLike TS>
void SumOperator<TS>::product(
    typename TS::Value inScale,
    typename TS::ConstVectorView in,
    typename TS::Value outScale,
    typename TS::VectorView out) {
        if (inScale) {
            _pLhs->product(0, in, 1, _buffer);
            _pRhs->product(1, in, 1, _buffer);
            _pSpace->scale(out, inScale);
            _pSpace->add(out, _buffer, outScale);
        } else {
            _pLhs->product(0, in, 1, out);
            _pRhs->product(outScale, in, outScale, out);
        }
}


} // namespace cie::linalg
