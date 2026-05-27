#pragma once

// --- External Includes ---
#include <Eigen/Dense> // Eigen::Map

// help the language server
#include "packages/integrands/inc/DirichletPenaltyIntegrand.hpp"


namespace cie::fem {


template <
    maths::Expression TDirichlet,
    maths::Expression TAnsatz,
    maths::Expression TEmbedding,
    CellLike TCell,
    bool Symmetric>
DirichletPenaltyIntegrand<TDirichlet,TAnsatz,TEmbedding,TCell,Symmetric>::DirichletPenaltyIntegrand() noexcept
    : _penalty(0),
      _dirichletFunctor(),
      _ansatzSpace(),
      _embedding(),
      _cellInverseTransform(),
      _cellJacobian()
{}


template <
    maths::Expression TDirichlet,
    maths::Expression TAnsatz,
    maths::Expression TEmbedding,
    CellLike TCell,
    bool Symmetric>
DirichletPenaltyIntegrand<TDirichlet,TAnsatz,TEmbedding,TCell,Symmetric>::DirichletPenaltyIntegrand(
    Ref<const TDirichlet> rDirichletFunctor,
    const Value penalty,
    Ref<const TAnsatz> rAnsatzSpace,
    Ref<const TEmbedding> rBoundaryTransform,
    Ref<const CellInverseTransform> rCellInverseTransform,
    Ref<const CellJacobian> rCellJacobian)
        : _penalty(penalty),
          _dirichletFunctor(rDirichletFunctor),
          _ansatzSpace(rAnsatzSpace),
          _embedding(rBoundaryTransform),
          _cellInverseTransform(rCellInverseTransform),
          _cellJacobian(rCellJacobian)
{}


template <
    maths::Expression TDirichlet,
    maths::Expression TAnsatz,
    maths::Expression TEmbedding,
    CellLike TCell,
    bool Symmetric>
void DirichletPenaltyIntegrand<TDirichlet,TAnsatz,TEmbedding,TCell,Symmetric>::evaluate(
    ConstSpan in,
    Span out,
    BufferSpan buffer) const {
        assert(this->bufferSize() <= buffer.size());

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

        const Span
            ansatzBuffer(pAnsatzBufferBegin, ansatzCount),
            boundaryTransformedBuffer(pEmbeddingBufferBegin, TCell::PhysicalDimension),
            cellTransformedBuffer(pCellParametricBufferBegin, TCell::ParametricDimension),
            dirichletBuffer(pDirichletBufferBegin, stateVariableCount),
            nestedBuffer(pNestedBufferBegin, buffer.data() + buffer.size());

        // Compute the cell's jacobian's determinant at the given input.
        const Value jacobian = std::abs(_cellJacobian.evaluateDeterminant(
            in,
            {}));
        const Value scaledPenalty = jacobian * _penalty;

        // Compute LHS contribution.
        _embedding.evaluate(in, boundaryTransformedBuffer, nestedBuffer);
        _cellInverseTransform.evaluate(boundaryTransformedBuffer, cellTransformedBuffer, nestedBuffer);

        if constexpr (Symmetric) {
            _ansatzSpace.evaluate(cellTransformedBuffer, ansatzBuffer, nestedBuffer);
            EigenAdaptor ansatzAdaptor(ansatzBuffer.data(), ansatzBuffer.size(), 1);
            EigenAdaptor lhsAdaptor(out.data(), ansatzCount, ansatzCount);
            lhsAdaptor.noalias() = ansatzAdaptor * scaledPenalty * ansatzAdaptor.transpose();
        } /*Symmetric*/ else {
            _ansatzSpace.evaluate(
                cellTransformedBuffer,
                {out.data(), ansatzCount},
                nestedBuffer);
        }

        // Compute the prescribed state at the given input.
        _dirichletFunctor.evaluate(in, dirichletBuffer, nestedBuffer);

        // Compute RHS contribution.
        for (unsigned iStateComponent=0u; iStateComponent<dirichletBuffer.size(); ++iStateComponent) {
            const Value dirichlet = dirichletBuffer[iStateComponent];
            const std::size_t rhsOffset = Symmetric
                ? ansatzCount * ansatzCount
                : ansatzCount;
            std::transform(
                ansatzBuffer.begin(),
                ansatzBuffer.end(),
                out.begin() + rhsOffset + iStateComponent * ansatzCount,
                [scaledPenalty, dirichlet] (const Value ansatz) -> Value {
                    return dirichlet * ansatz * scaledPenalty;
                });
        } // for iStateComponent
}


template <
    maths::Expression TDirichlet,
    maths::Expression TAnsatz,
    maths::Expression TEmbedding,
    CellLike TCell,
    bool Symmetric>
unsigned DirichletPenaltyIntegrand<TDirichlet,TAnsatz,TEmbedding,TCell,Symmetric>::size() const {
    const auto ansatzCount = _ansatzSpace.size();
    return
          ansatzCount * ansatzCount     // <= LHS contribution
        + ansatzCount;                  // <= RHS contribution
}


template <
    maths::Expression TDirichlet,
    maths::Expression TAnsatz,
    maths::Expression TEmbedding,
    CellLike TCell,
    bool Symmetric>
unsigned DirichletPenaltyIntegrand<TDirichlet,TAnsatz,TEmbedding,TCell,Symmetric>::bufferSize() const noexcept {
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
