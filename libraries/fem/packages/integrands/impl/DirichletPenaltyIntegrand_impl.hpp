#pragma once

// --- External Includes ---
#include <Eigen/Dense> // Eigen::Map

// --- FEM Includes ---
#include "packages/integrands/inc/DirichletPenaltyIntegrand.hpp"

// --- Utility Includes ---
#include "packages/io/inc/Serializer.hpp"


namespace cie::fem {


template <
    maths::Expression TDirichlet,
    maths::Expression TAnsatz,
    maths::Expression TEmbedding,
    CellLike TCell>
DirichletPenaltyIntegrand<TDirichlet,TAnsatz,TEmbedding,TCell>::DirichletPenaltyIntegrand() noexcept
    :   _ansatzSpace(),
        _embedding(),
        _cellInverseTransform(),
        _cellJacobian(),
        _dirichletFunctor(),
        _penalty(0)
{}


template <
    maths::Expression TDirichlet,
    maths::Expression TAnsatz,
    maths::Expression TEmbedding,
    CellLike TCell>
DirichletPenaltyIntegrand<TDirichlet,TAnsatz,TEmbedding,TCell>::DirichletPenaltyIntegrand(
    Ref<const TDirichlet> rDirichletFunctor,
    const Value penalty,
    Ref<const TAnsatz> rAnsatzSpace,
    Ref<const TEmbedding> rBoundaryTransform,
    Ref<const CellInverseTransform> rCellInverseTransform,
    Ref<const CellJacobian> rCellJacobian)
        :   _ansatzSpace(rAnsatzSpace),
            _embedding(rBoundaryTransform),
            _cellInverseTransform(rCellInverseTransform),
            _cellJacobian(rCellJacobian),
            _dirichletFunctor(rDirichletFunctor),
            _penalty(penalty)
{}


template <
    maths::Expression TDirichlet,
    maths::Expression TAnsatz,
    maths::Expression TEmbedding,
    CellLike TCell>
void DirichletPenaltyIntegrand<TDirichlet,TAnsatz,TEmbedding,TCell>::evaluate(
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

        _ansatzSpace.evaluate(cellTransformedBuffer, ansatzBuffer, nestedBuffer);
        EigenAdaptor ansatzAdaptor(ansatzBuffer.data(), ansatzBuffer.size(), 1);
        EigenAdaptor lhsAdaptor(out.data(), ansatzCount, ansatzCount);
        lhsAdaptor.noalias() = ansatzAdaptor * scaledPenalty * ansatzAdaptor.transpose();

        // Compute the prescribed state at the given input.
        _dirichletFunctor.evaluate(in, dirichletBuffer, nestedBuffer);

        // Compute RHS contribution.
        for (unsigned iStateComponent=0u; iStateComponent<dirichletBuffer.size(); ++iStateComponent) {
            const Value dirichlet = dirichletBuffer[iStateComponent];
            const std::size_t rhsOffset = ansatzCount * ansatzCount;
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
    CellLike TCell>
unsigned DirichletPenaltyIntegrand<TDirichlet,TAnsatz,TEmbedding,TCell>::size() const {
    const auto ansatzCount = _ansatzSpace.size();
    return
          ansatzCount * ansatzCount     // <= LHS contribution
        + ansatzCount;                  // <= RHS contribution
}


template <
    maths::Expression TDirichlet,
    maths::Expression TAnsatz,
    maths::Expression TEmbedding,
    CellLike TCell>
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


template <
    maths::Expression TDirichlet,
    maths::Expression TAnsatz,
    maths::Expression TEmbedding,
    CellLike TCell>
void DirichletPenaltyIntegrand<TDirichlet,TAnsatz,TEmbedding,TCell>::serialize(
    Ref<cie::io::Traits::SerializerStream> rStream,
    tags::Binary) const {
        using BS = cie::io::BinarySerializer;
        BS::serialize(
            rStream,
            _ansatzSpace);
        BS::serialize(
            rStream,
            _embedding);
        BS::serialize(
            rStream,
            _cellInverseTransform);
        BS::serialize(
            rStream,
            _cellJacobian);
        BS::serialize(
            rStream,
            _dirichletFunctor);
        BS::serialize(
            rStream,
            _penalty);
}


template <
    maths::Expression TDirichlet,
    maths::Expression TAnsatz,
    maths::Expression TEmbedding,
    CellLike TCell>
template <class TAllocator>
void DirichletPenaltyIntegrand<TDirichlet,TAnsatz,TEmbedding,TCell>::deserialize(
    Ref<cie::io::Traits::DeserializerStream> rStream,
    Ref<DirichletPenaltyIntegrand> rInstance,
    Ref<const TAllocator> rAllocator,
    tags::Binary) {
        using BS = cie::io::BinarySerializer;
        BS::deserialize(
            rStream,
            rInstance._ansatzSpace,
            rAllocator);
        BS::deserialize(
            rStream,
            rInstance._embedding,
            rAllocator);
        BS::deserialize(
            rStream,
            rInstance._cellInverseTransform,
            rAllocator);
        BS::deserialize(
            rStream,
            rInstance._cellJacobian,
            rAllocator);
        BS::deserialize(
            rStream,
            rInstance._dirichletFunctor,
            rAllocator);
        BS::deserialize(
            rStream,
            rInstance._penalty,
            rAllocator);
}


} // namespace cie::fem
