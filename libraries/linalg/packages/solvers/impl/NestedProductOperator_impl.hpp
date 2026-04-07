#pragma once

// --- Linalg Includes ---
#include "packages/solvers/inc/NestedProductOperator.hpp"


namespace cie::linalg {


template <LinalgSpaceLike TS>
NestedProductOperator<TS>::NestedProductOperator(
    std::shared_ptr<TS> pSpace,
    std::shared_ptr<LinearOperator<TS>> pLhs,
    std::shared_ptr<LinearOperator<TS>> pRhs,
    std::size_t rowCount)
        :   _pSpace(pSpace),
            _pLhs(pLhs),
            _pRhs(pRhs),
            _buffer(pSpace->makeVector(rowCount))
{}


template <LinalgSpaceLike TS>
void NestedProductOperator<TS>::product(
    typename TS::Value inScale,
    typename TS::ConstVectorView in,
    typename TS::Value outScale,
    typename TS::VectorView out) {
        _pRhs->product(0, in, 1, _buffer);
        _pLhs->product(inScale, _buffer, outScale, out);
}


} // namespace cie::linalg
