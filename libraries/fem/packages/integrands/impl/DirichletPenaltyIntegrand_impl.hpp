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
      _cellInverseTransform()
{}


template <maths::Expression TDirichlet, maths::Expression TAnsatz, maths::Expression TEmbedding, CellLike TCell>
DirichletPenaltyIntegrand<TDirichlet,TAnsatz,TEmbedding,TCell>::DirichletPenaltyIntegrand(
    Ref<const TDirichlet> rDirichletFunctor,
    const Value penalty,
    Ref<const TAnsatz> rAnsatzSpace,
    Ref<const TEmbedding> rBoundaryTransform,
    Ref<const CellInverseTransform> rCellInverseTransform)
        : _penalty(penalty),
          _dirichletFunctor(rDirichletFunctor),
          _ansatzSpace(rAnsatzSpace),
          _embedding(rBoundaryTransform),
          _cellInverseTransform(rCellInverseTransform)
{}


template <maths::Expression TDirichlet, maths::Expression TAnsatz, maths::Expression TEmbedding, CellLike TCell>
void DirichletPenaltyIntegrand<TDirichlet,TAnsatz,TEmbedding,TCell>::evaluate(
    ConstSpan in,
    Span out,
    BufferSpan buffer) const {
        assert(this->getMinBufferSize() <= buffer.size());

        // Fetch object output sizes.
        const unsigned ansatzCount = _ansatzSpace.size();
        const unsigned stateVariableCount = _dirichletFunctor.size();

        using EigenDenseMatrix = Eigen::Matrix<Value,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>;
        using EigenAdaptor = Eigen::Map<EigenDenseMatrix>;

        Ptr<Value> pAnsatzBufferBegin = buffer.data();
        Ptr<Value> pEmbeddingBufferBegin = pAnsatzBufferBegin + ansatzCount;
        Ptr<Value> pCellParametricBufferBegin = pEmbeddingBufferBegin + TCell::PhysicalDimension;
        Ptr<Value> pDirichletBufferBegin = pCellParametricBufferBegin + TCell::ParametricDimension;
        Ptr<Value> pNestedBufferBegin = pDirichletBufferBegin + stateVariableCount;

        const Span ansatzBuffer(pAnsatzBufferBegin, ansatzCount),
                boundaryTransformedBuffer(pEmbeddingBufferBegin, TCell::PhysicalDimension),
                cellTransformedBuffer(pCellParametricBufferBegin, TCell::ParametricDimension),
                dirichletBuffer(pDirichletBufferBegin, stateVariableCount),
                nestedBuffer(pNestedBufferBegin, buffer.data() + buffer.size());

        // Compute LHS contribution.
        _embedding.evaluate(in, boundaryTransformedBuffer, nestedBuffer);
        _cellInverseTransform.evaluate(boundaryTransformedBuffer, cellTransformedBuffer, nestedBuffer);

        _ansatzSpace.evaluate(cellTransformedBuffer, ansatzBuffer, nestedBuffer);
        EigenAdaptor ansatzAdaptor(ansatzBuffer.data(), ansatzBuffer.size(), 1);
        EigenAdaptor lhsAdaptor(out.data(), ansatzCount, ansatzCount);

        lhsAdaptor.noalias() = ansatzAdaptor * _penalty * ansatzAdaptor.transpose();

        // Compute the prescribed state at the given input.
        _dirichletFunctor.evaluate(in, dirichletBuffer, nestedBuffer);

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
unsigned DirichletPenaltyIntegrand<TDirichlet,TAnsatz,TEmbedding,TCell>::bufferSize() const noexcept {
    const std::array<unsigned,4> nestedBufferRequirements {
        _dirichletFunctor.bufferSize(),
        _ansatzSpace.bufferSize(),
        _embedding.bufferSize(),
        _cellInverseTransform.bufferSize()};
    const unsigned nestedBufferSize = *std::max_element(
        nestedBufferRequirements.begin(),
        nestedBufferRequirements.end());
    return _ansatzSpace.size()
         + TCell::PhysicalDimension
         + TCell::ParametricDimension
         + _dirichletFunctor.size()
         + nestedBufferSize;
}


} // namespace cie::fem
