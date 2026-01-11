#pragma once

// --- External Includes ---
#include <Eigen/Dense> // Eigen::Map

// help the language server
#include "packages/integrands/inc/DirichletPenaltyIntegrand.hpp"


namespace cie::fem {


template <maths::Expression TDirichlet, maths::Expression TAnsatz, maths::Expression TTransform>
DirichletPenaltyIntegrand<TDirichlet,TAnsatz,TTransform>::DirichletPenaltyIntegrand()
    : DirichletPenaltyIntegrand(TDirichlet(), 0, nullptr)
{
}


template <maths::Expression TDirichlet, maths::Expression TAnsatz, maths::Expression TTransform>
DirichletPenaltyIntegrand<TDirichlet,TAnsatz,TTransform>::DirichletPenaltyIntegrand(Ref<const TDirichlet> rDirichletFunctor,
                                                                                    const Value penalty,
                                                                                    Ref<const TAnsatz> rAnsatzSpace,
                                                                                    Ref<const TTransform> rSpatialTransform)
    : _penalty(penalty),
      _pDirichletFunctor(&rDirichletFunctor),
      _pAnsatzSpace(&rAnsatzSpace),
      _pSpatialTransform(&rSpatialTransform),
      _buffer()
{
}


template <maths::Expression TDirichlet, maths::Expression TAnsatz, maths::Expression TTransform>
DirichletPenaltyIntegrand<TDirichlet,TAnsatz,TTransform>::DirichletPenaltyIntegrand(Ref<const TDirichlet> rDirichletFunctor,
                                                                                    const Value penalty,
                                                                                    Ref<const TAnsatz> rAnsatzSpace,
                                                                                    Ref<const TTransform> rSpatialTransform,
                                                                                    std::span<Value> buffer)
    : DirichletPenaltyIntegrand(rDirichletFunctor, penalty, rAnsatzSpace, rSpatialTransform)
{
    this->setBuffer(buffer);
}


template <maths::Expression TDirichlet, maths::Expression TAnsatz, maths::Expression TTransform>
void DirichletPenaltyIntegrand<TDirichlet,TAnsatz,TTransform>::evaluate(ConstSpan in, Span out) const {
    CIE_OUT_OF_RANGE_CHECK(this->getMinBufferSize() <= _buffer.size())

    Ref<const TAnsatz> rAnsatzSpace = *_pAnsatzSpace;

    // Fetch object output sizes.
    const unsigned ansatzCount = rAnsatzSpace.size();
    const unsigned transformedSize = _pSpatialTransform->size();
    const unsigned stateVariableCount = _pDirichletFunctor->size();

    using EigenDenseMatrix = Eigen::Matrix<Value,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>;
    using EigenAdaptor = Eigen::Map<EigenDenseMatrix>;

    auto pAnsatzBufferBegin = _buffer.data();
    auto pTransformedBufferBegin = pAnsatzBufferBegin + ansatzCount;
    auto pDirichletBufferBegin = pTransformedBufferBegin + transformedSize;

    Span ansatzBuffer(pAnsatzBufferBegin, ansatzCount),
         transformedBuffer(pTransformedBufferBegin, transformedSize),
         dirichletBuffer(pDirichletBufferBegin, stateVariableCount);

    // Compute LHS contribution.
    rAnsatzSpace.evaluate(in, ansatzBuffer);
    EigenAdaptor ansatzAdaptor(ansatzBuffer.data(), ansatzBuffer.size(), 1);
    EigenAdaptor lhsAdaptor(out.data(), ansatzCount, ansatzCount);

    lhsAdaptor = ansatzAdaptor * _penalty * ansatzAdaptor.transpose();

    // Compute the prescribed state at the given input.
    _pSpatialTransform->evaluate(in, transformedBuffer);
    _pDirichletFunctor->evaluate(transformedBuffer, dirichletBuffer);

    // Compute RHS contribution.
    for (unsigned iStateComponent=0u; iStateComponent<dirichletBuffer.size(); ++iStateComponent) {
        const Value dirichlet = dirichletBuffer[iStateComponent];
        std::transform(
            ansatzBuffer.begin(),
            ansatzBuffer.end(),
            out.begin() + ansatzCount * ansatzCount + iStateComponent * ansatzCount,
            [this, dirichlet] (const Value ansatz) -> Value {
                return dirichlet * ansatz * _penalty;
            });
    }
}


template <maths::Expression TDirichlet, maths::Expression TAnsatz, maths::Expression TTransform>
unsigned DirichletPenaltyIntegrand<TDirichlet,TAnsatz,TTransform>::size() const {
    const auto ansatzCount = _pAnsatzSpace->size();
    return
          ansatzCount * ansatzCount     // <= LHS contribution
        + ansatzCount;                  // <= RHS contribution
}


template <maths::Expression TDirichlet, maths::Expression TAnsatz, maths::Expression TTransform>
void DirichletPenaltyIntegrand<TDirichlet,TAnsatz,TTransform>::setBuffer(std::span<Value> buffer) {
    CIE_OUT_OF_RANGE_CHECK(this->getMinBufferSize() <= buffer.size())
    _buffer = buffer;
}


template <maths::Expression TDirichlet, maths::Expression TAnsatz, maths::Expression TTransform>
unsigned DirichletPenaltyIntegrand<TDirichlet,TAnsatz,TTransform>::getMinBufferSize() const noexcept {
    return _pAnsatzSpace->size() + _pSpatialTransform->size() + _pDirichletFunctor->size();
}


} // namespace cie::fem
