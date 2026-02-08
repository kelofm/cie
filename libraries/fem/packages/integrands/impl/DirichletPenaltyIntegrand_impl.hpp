#pragma once

// --- External Includes ---
#include <Eigen/Dense> // Eigen::Map

// help the language server
#include "packages/integrands/inc/DirichletPenaltyIntegrand.hpp"


namespace cie::fem {


template <maths::Expression TDirichlet, maths::Expression TAnsatz, maths::Expression TEmbedding, CellLike TCell>
DirichletPenaltyIntegrand<TDirichlet,TAnsatz,TEmbedding,TCell>::DirichletPenaltyIntegrand()
    : _penalty(0),
      _dirichletFunctor(),
      _ansatzSpace(),
      _embedding(),
      _cellInverseTransform(),
      _buffer()
{}


template <maths::Expression TDirichlet, maths::Expression TAnsatz, maths::Expression TEmbedding, CellLike TCell>
DirichletPenaltyIntegrand<TDirichlet,TAnsatz,TEmbedding,TCell>::DirichletPenaltyIntegrand(Ref<const TDirichlet> rDirichletFunctor,
                                                                                         const Value penalty,
                                                                                         Ref<const TAnsatz> rAnsatzSpace,
                                                                                         Ref<const TEmbedding> rBoundaryTransform,
                                                                                         Ref<const CellInverseTransform> rCellInverseTransform)
    : _penalty(penalty),
      _dirichletFunctor(rDirichletFunctor),
      _ansatzSpace(rAnsatzSpace),
      _embedding(rBoundaryTransform),
      _cellInverseTransform(rCellInverseTransform),
      _buffer()
{}


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

    // Fetch object output sizes.
    const unsigned ansatzCount = _ansatzSpace.size();
    const unsigned stateVariableCount = _dirichletFunctor.size();

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
    _embedding.evaluate(in, boundaryTransformedBuffer);
    _cellInverseTransform.evaluate(boundaryTransformedBuffer, cellTransformedBuffer);

    _ansatzSpace.evaluate(cellTransformedBuffer, ansatzBuffer);
    EigenAdaptor ansatzAdaptor(ansatzBuffer.data(), ansatzBuffer.size(), 1);
    EigenAdaptor lhsAdaptor(out.data(), ansatzCount, ansatzCount);

    lhsAdaptor.noalias() = ansatzAdaptor * _penalty * ansatzAdaptor.transpose();

    // Compute the prescribed state at the given input.
    _dirichletFunctor.evaluate(boundaryTransformedBuffer, dirichletBuffer);

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
    const auto ansatzCount = _ansatzSpace.size();
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
    return _ansatzSpace.size()
         + TCell::PhysicalDimension
         + TCell::ParametricDimension
         + _dirichletFunctor.size();
}


} // namespace cie::fem
