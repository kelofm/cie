#pragma once

// --- External Includes ---
#include <Eigen/Dense> // Eigen::Map

// help the language server
#include "packages/integrands/inc/DirichletPenaltyIntegrand.hpp"


namespace cie::fem {


template <maths::Expression TDirichlet, maths::Expression TAnsatz, maths::Expression TEmbedding, CellLike TCell>
DirichletPenaltyIntegrand<TDirichlet,TAnsatz,TEmbedding,TCell>::DirichletPenaltyIntegrand()
    : DirichletPenaltyIntegrand(TDirichlet(), 0, nullptr)
{
}


template <maths::Expression TDirichlet, maths::Expression TAnsatz, maths::Expression TEmbedding, CellLike TCell>
DirichletPenaltyIntegrand<TDirichlet,TAnsatz,TEmbedding,TCell>::DirichletPenaltyIntegrand(Ref<const TDirichlet> rDirichletFunctor,
                                                                                         const Value penalty,
                                                                                         Ref<const TAnsatz> rAnsatzSpace,
                                                                                         Ref<const TEmbedding> rBoundaryTransform,
                                                                                         Ref<const CellInverseTransform> rCellInverseTransform)
    : _penalty(penalty),
      _pDirichletFunctor(&rDirichletFunctor),
      _pAnsatzSpace(&rAnsatzSpace),
      _pEmbedding(&rBoundaryTransform),
      _pCellInverseTransform(&rCellInverseTransform),
      _buffer()
{
}


template <maths::Expression TDirichlet, maths::Expression TAnsatz, maths::Expression TEmbedding, CellLike TCell>
DirichletPenaltyIntegrand<TDirichlet,TAnsatz,TEmbedding,TCell>::DirichletPenaltyIntegrand(Ref<const TDirichlet> rDirichletFunctor,
                                                                                         const Value penalty,
                                                                                         Ref<const TAnsatz> rAnsatzSpace,
                                                                                         Ref<const TEmbedding> rEmbedding,
                                                                                         Ref<const CellInverseTransform> rCellInverseTransform,
                                                                                         std::span<Value> buffer)
    : DirichletPenaltyIntegrand(rDirichletFunctor, penalty, rAnsatzSpace, rEmbedding, rCellInverseTransform)
{
    this->setBuffer(buffer);
}


template <maths::Expression TDirichlet, maths::Expression TAnsatz, maths::Expression TEmbedding, CellLike TCell>
void DirichletPenaltyIntegrand<TDirichlet,TAnsatz,TEmbedding,TCell>::evaluate(ConstSpan in, Span out) const {
    CIE_OUT_OF_RANGE_CHECK(this->getMinBufferSize() <= _buffer.size())

    Ref<const TAnsatz> rAnsatzSpace = *_pAnsatzSpace;

    // Fetch object output sizes.
    const unsigned ansatzCount = rAnsatzSpace.size();
    const unsigned stateVariableCount = _pDirichletFunctor->size();

    using EigenDenseMatrix = Eigen::Matrix<Value,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>;
    using EigenAdaptor = Eigen::Map<EigenDenseMatrix>;

    auto pAnsatzBufferBegin = _buffer.data();
    auto pEmbeddingBufferBegin = pAnsatzBufferBegin + ansatzCount;
    auto pCellParametricBufferBegin = pEmbeddingBufferBegin + TCell::PhysicalDimension;
    auto pDirichletBufferBegin = pCellParametricBufferBegin + TCell::ParametricDimension;

    const Span ansatzBuffer(pAnsatzBufferBegin, ansatzCount),
               boundaryTransformedBuffer(pEmbeddingBufferBegin, TCell::PhysicalDimension),
               cellTransformedBuffer(pCellParametricBufferBegin, TCell::ParametricDimension),
               dirichletBuffer(pDirichletBufferBegin, stateVariableCount);

    // Compute LHS contribution.
    _pEmbedding->evaluate(in, boundaryTransformedBuffer);
    _pCellInverseTransform->evaluate(boundaryTransformedBuffer, cellTransformedBuffer);

    rAnsatzSpace.evaluate(cellTransformedBuffer, ansatzBuffer);
    EigenAdaptor ansatzAdaptor(ansatzBuffer.data(), ansatzBuffer.size(), 1);
    EigenAdaptor lhsAdaptor(out.data(), ansatzCount, ansatzCount);

    lhsAdaptor.noalias() = ansatzAdaptor * _penalty * ansatzAdaptor.transpose();

    // Compute the prescribed state at the given input.
    _pDirichletFunctor->evaluate(boundaryTransformedBuffer, dirichletBuffer);

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


template <maths::Expression TDirichlet, maths::Expression TAnsatz, maths::Expression TEmbedding, CellLike TCell>
unsigned DirichletPenaltyIntegrand<TDirichlet,TAnsatz,TEmbedding,TCell>::size() const {
    const auto ansatzCount = _pAnsatzSpace->size();
    return
          ansatzCount * ansatzCount     // <= LHS contribution
        + ansatzCount;                  // <= RHS contribution
}


template <maths::Expression TDirichlet, maths::Expression TAnsatz, maths::Expression TEmbedding, CellLike TCell>
void DirichletPenaltyIntegrand<TDirichlet,TAnsatz,TEmbedding,TCell>::setBuffer(std::span<Value> buffer) {
    CIE_OUT_OF_RANGE_CHECK(this->getMinBufferSize() <= buffer.size())
    _buffer = buffer;
}


template <maths::Expression TDirichlet, maths::Expression TAnsatz, maths::Expression TEmbedding, CellLike TCell>
unsigned DirichletPenaltyIntegrand<TDirichlet,TAnsatz,TEmbedding,TCell>::getMinBufferSize() const noexcept {
    return _pAnsatzSpace->size()
         + TCell::PhysicalDimension
         + TCell::ParametricDimension
         + _pDirichletFunctor->size();
}


} // namespace cie::fem
