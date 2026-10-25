#pragma once

// --- External Includes ---
#include <Eigen/Dense> // Eigen::Map
#include <algorithm>

// --- FEM Includes ---
#include "packages/integrands/inc/DirichletPenaltyIntegrand.hpp"

// --- Utility Includes ---
#include "packages/io/inc/Serializer.hpp"


namespace cie::fem {


template <
    maths::Expression TDirichlet,
    maths::SliceableExpression TAnsatz,
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
    maths::SliceableExpression TAnsatz,
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


/// @name Runtime Expression Interface
/// @{


template <
    maths::Expression TDirichlet,
    maths::SliceableExpression TAnsatz,
    maths::Expression TEmbedding,
    CellLike TCell>
void DirichletPenaltyIntegrand<TDirichlet,TAnsatz,TEmbedding,TCell>::evaluate(
    ConstSpan in,
    Span out,
    BufferSpan buffer) const {
        const auto bounds = this->contextBounds();
        Context context;
        for (unsigned iDimension=0u; iDimension<ParametricDimension; ++iDimension) {
            context.data()[2 * iDimension] = 0;
            context.data()[2 * iDimension + 1] = bounds[iDimension];
        }
        return this->evaluate(
            in,
            out,
            buffer,
            context);
}


template <
    maths::Expression TDirichlet,
    maths::SliceableExpression TAnsatz,
    maths::Expression TEmbedding,
    CellLike TCell>
void DirichletPenaltyIntegrand<TDirichlet,TAnsatz,TEmbedding,TCell>::evaluate(
    ConstSpan in,
    Span out,
    BufferSpan buffer,
    Ref<const Context> rContext) const {
        assert(this->bufferSize(rContext) <= buffer.size());

        // Fetch object output sizes.
        const unsigned ansatzCount = _ansatzSpace.size(rContext);
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

        _ansatzSpace.evaluate(cellTransformedBuffer, ansatzBuffer, nestedBuffer, rContext);
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
    maths::SliceableExpression TAnsatz,
    maths::Expression TEmbedding,
    CellLike TCell>
unsigned DirichletPenaltyIntegrand<TDirichlet,TAnsatz,TEmbedding,TCell>::size() const noexcept {
    const auto ansatzCount = _ansatzSpace.size();
    return
          ansatzCount * ansatzCount     // <= LHS contribution
        + ansatzCount;                  // <= RHS contribution
}


template <
    maths::Expression TDirichlet,
    maths::SliceableExpression TAnsatz,
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
    maths::SliceableExpression TAnsatz,
    maths::Expression TEmbedding,
    CellLike TCell>
unsigned DirichletPenaltyIntegrand<TDirichlet,TAnsatz,TEmbedding,TCell>::size(Ref<const Context> rContext) const noexcept {
    const auto ansatzCount = _ansatzSpace.size(rContext);
    return
          ansatzCount * ansatzCount     // <= LHS contribution
        + ansatzCount;                  // <= RHS contribution
}


template <
    maths::Expression TDirichlet,
    maths::SliceableExpression TAnsatz,
    maths::Expression TEmbedding,
    CellLike TCell>
unsigned DirichletPenaltyIntegrand<TDirichlet,TAnsatz,TEmbedding,TCell>::bufferSize(Ref<const Context> rContext) const noexcept {
    const std::array<unsigned,4> nestedBufferRequirements {
        _dirichletFunctor.bufferSize(),
        _ansatzSpace.bufferSize(rContext),
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
    maths::Expression TD,
    maths::SliceableExpression TA,
    maths::Expression TE,
    CellLike TC>
std::array<
    typename DirichletPenaltyIntegrand<TD,TA,TE,TC>::Context::Value,
    DirichletPenaltyIntegrand<TD,TA,TE,TC>::ParametricDimension>
DirichletPenaltyIntegrand<TD,TA,TE,TC>::contextBounds() const noexcept {
    return _ansatzSpace.contextBounds();
}


/// @}


template <
    maths::Expression TDirichlet,
    maths::SliceableExpression TAnsatz,
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
    maths::SliceableExpression TAnsatz,
    maths::Expression TEmbedding,
    CellLike TCell>
template <class TStream, class TAllocator>
void DirichletPenaltyIntegrand<TDirichlet,TAnsatz,TEmbedding,TCell>::deserialize(
    Ref<TStream> rStream,
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


template <
    maths::SliceableExpression TAnsatz,
    maths::Expression TEmbedding,
    CellLike TCell>
DirichletPenaltyLHSIntegrand<TAnsatz,TEmbedding,TCell>::DirichletPenaltyLHSIntegrand() noexcept
    :   _ansatzSpace(),
        _embedding(),
        _cellInverseTransform(),
        _cellJacobian(),
        _penalty(0)
{}


template <
    maths::SliceableExpression TAnsatz,
    maths::Expression TEmbedding,
    CellLike TCell>
DirichletPenaltyLHSIntegrand<TAnsatz,TEmbedding,TCell>::DirichletPenaltyLHSIntegrand(
    const Value penalty,
    Ref<const TAnsatz> rAnsatzSpace,
    Ref<const TEmbedding> rBoundaryTransform,
    Ref<const CellInverseTransform> rCellInverseTransform,
    Ref<const CellJacobian> rCellJacobian)
        :   _ansatzSpace(rAnsatzSpace),
            _embedding(rBoundaryTransform),
            _cellInverseTransform(rCellInverseTransform),
            _cellJacobian(rCellJacobian),
            _penalty(penalty)
{}


template <
    maths::SliceableExpression TAnsatz,
    maths::Expression TEmbedding,
    CellLike TCell>
void DirichletPenaltyLHSIntegrand<TAnsatz,TEmbedding,TCell>::evaluate(
    ConstSpan in,
    Span out,
    BufferSpan buffer) const {
        Context context;
        const auto bounds = this->contextBounds();
        for (unsigned iDimension=0u; iDimension<ParametricDimension; ++iDimension) {
            context.data()[2 * iDimension] = 0;
            context.data()[2 * iDimension + 1] = bounds[iDimension];
        }
        const InContext inContext{context,context};
        return this->evaluate(in, out, buffer, inContext);
}


template <
    maths::SliceableExpression TAnsatz,
    maths::Expression TEmbedding,
    CellLike TCell>
void DirichletPenaltyLHSIntegrand<TAnsatz,TEmbedding,TCell>::evaluate(
    ConstSpan in,
    Span out,
    BufferSpan buffer,
    Ref<const InContext> rInContext) const {
        using AnsatzInContext = typename TAnsatz::InContext;
        const AnsatzInContext rowContext{rInContext[0]};
        const AnsatzInContext columnContext{rInContext[1]};
        const unsigned rowCount = _ansatzSpace.size(rowContext);
        const unsigned columnCount = _ansatzSpace.size(columnContext);
        assert(out.size() == rowCount * columnCount);
        assert(this->bufferSize(rInContext) <= buffer.size());

        Ptr<Value> pRowAnsatzBuffer = buffer.data();
        Ptr<Value> pColumnAnsatzBuffer = pRowAnsatzBuffer + rowCount;
        Ptr<Value> pEmbeddingBuffer = pColumnAnsatzBuffer + columnCount;
        Ptr<Value> pCellParametricBuffer = pEmbeddingBuffer + TCell::PhysicalDimension;
        Ptr<Value> pNestedBuffer = pCellParametricBuffer + TCell::ParametricDimension;
        const Span rowAnsatzBuffer(pRowAnsatzBuffer, rowCount);
        const Span columnAnsatzBuffer(pColumnAnsatzBuffer, columnCount);
        const Span boundaryTransformedBuffer(pEmbeddingBuffer, TCell::PhysicalDimension);
        const Span cellTransformedBuffer(pCellParametricBuffer, TCell::ParametricDimension);
        const BufferSpan nestedBuffer(pNestedBuffer, buffer.data() + buffer.size());

        const Value jacobian = std::abs(_cellJacobian.evaluateDeterminant(in, {}));
        const Value scaledPenalty = jacobian * _penalty;
        _embedding.evaluate(in, boundaryTransformedBuffer, nestedBuffer);
        _cellInverseTransform.evaluate(boundaryTransformedBuffer, cellTransformedBuffer, nestedBuffer);
        _ansatzSpace.evaluate(cellTransformedBuffer, rowAnsatzBuffer, nestedBuffer, rowContext);
        _ansatzSpace.evaluate(cellTransformedBuffer, columnAnsatzBuffer, nestedBuffer, columnContext);

        using DenseMatrix = Eigen::Matrix<Value,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>;
        using MatrixAdaptor = Eigen::Map<DenseMatrix>;
        MatrixAdaptor rowAdaptor(rowAnsatzBuffer.data(), rowCount, 1);
        MatrixAdaptor columnAdaptor(columnAnsatzBuffer.data(), columnCount, 1);
        MatrixAdaptor lhsAdaptor(out.data(), rowCount, columnCount);
        lhsAdaptor.noalias() = rowAdaptor * scaledPenalty * columnAdaptor.transpose();
}


template <
    maths::SliceableExpression TAnsatz,
    maths::Expression TEmbedding,
    CellLike TCell>
unsigned DirichletPenaltyLHSIntegrand<TAnsatz,TEmbedding,TCell>::size() const noexcept {
    const unsigned ansatzCount = _ansatzSpace.size();
    return ansatzCount * ansatzCount;
}


template <
    maths::SliceableExpression TAnsatz,
    maths::Expression TEmbedding,
    CellLike TCell>
unsigned DirichletPenaltyLHSIntegrand<TAnsatz,TEmbedding,TCell>::size(
    Ref<const InContext> rInContext) const noexcept {
        using AnsatzInContext = typename TAnsatz::InContext;
        return _ansatzSpace.size(AnsatzInContext{rInContext[0]})
             * _ansatzSpace.size(AnsatzInContext{rInContext[1]});
}


template <
    maths::SliceableExpression TAnsatz,
    maths::Expression TEmbedding,
    CellLike TCell>
unsigned DirichletPenaltyLHSIntegrand<TAnsatz,TEmbedding,TCell>::bufferSize() const noexcept {
    const auto bounds = this->contextBounds();
    Context context;
    for (unsigned iDimension=0u; iDimension<ParametricDimension; ++iDimension) {
        context.data()[2 * iDimension] = 0;
        context.data()[2 * iDimension + 1] = bounds[iDimension];
    }
    const InContext inContext{context,context};
    return this->bufferSize(inContext);
}


template <
    maths::SliceableExpression TAnsatz,
    maths::Expression TEmbedding,
    CellLike TCell>
unsigned DirichletPenaltyLHSIntegrand<TAnsatz,TEmbedding,TCell>::bufferSize(
    Ref<const InContext> rInContext) const noexcept {
        using AnsatzInContext = typename TAnsatz::InContext;
        const AnsatzInContext rowContext{rInContext[0]};
        const AnsatzInContext columnContext{rInContext[1]};
        const unsigned rowCount = _ansatzSpace.size(rowContext);
        const unsigned columnCount = _ansatzSpace.size(columnContext);
        const std::array<unsigned,4> nestedBufferRequirements {
            _ansatzSpace.bufferSize(rowContext),
            _ansatzSpace.bufferSize(columnContext),
            _embedding.bufferSize(),
            _cellInverseTransform.bufferSize()};
        return rowCount + columnCount
             + TCell::PhysicalDimension
             + TCell::ParametricDimension
             + *std::max_element(nestedBufferRequirements.begin(), nestedBufferRequirements.end());
}


template <
    maths::SliceableExpression TAnsatz,
    maths::Expression TEmbedding,
    CellLike TCell>
std::array<
    typename DirichletPenaltyLHSIntegrand<TAnsatz,TEmbedding,TCell>::Context::Value,
    DirichletPenaltyLHSIntegrand<TAnsatz,TEmbedding,TCell>::ParametricDimension>
DirichletPenaltyLHSIntegrand<TAnsatz,TEmbedding,TCell>::contextBounds() const noexcept {
    return _ansatzSpace.contextBounds();
}


template <
    maths::SliceableExpression TAnsatz,
    maths::Expression TEmbedding,
    CellLike TCell>
typename DirichletPenaltyLHSIntegrand<TAnsatz,TEmbedding,TCell>::OutContext
DirichletPenaltyLHSIntegrand<TAnsatz,TEmbedding,TCell>::outContext(
    Ref<const InContext> rInContext) const noexcept {
        std::array<typename OutContext::Value,4 * ParametricDimension> data{};
        for (unsigned iDimension=0u; iDimension<ParametricDimension; ++iDimension) {
            data[2 * iDimension] = rInContext[1].data()[2 * iDimension];
            data[2 * iDimension + 1] = rInContext[1].data()[2 * iDimension + 1];
            data[2 * (ParametricDimension + iDimension)] = rInContext[0].data()[2 * iDimension];
            data[2 * (ParametricDimension + iDimension) + 1] = rInContext[0].data()[2 * iDimension + 1];
        }
        return OutContext(data);
}


template <
    maths::SliceableExpression TAnsatz,
    maths::Expression TEmbedding,
    CellLike TCell>
void DirichletPenaltyLHSIntegrand<TAnsatz,TEmbedding,TCell>::serialize(
    Ref<cie::io::Traits::SerializerStream> rStream,
    tags::Binary) const {
        using BS = cie::io::BinarySerializer;
        BS::serialize(rStream, _ansatzSpace);
        BS::serialize(rStream, _embedding);
        BS::serialize(rStream, _cellInverseTransform);
        BS::serialize(rStream, _cellJacobian);
        BS::serialize(rStream, _penalty);
}


template <
    maths::SliceableExpression TAnsatz,
    maths::Expression TEmbedding,
    CellLike TCell>
template <class TStream, class TAllocator>
void DirichletPenaltyLHSIntegrand<TAnsatz,TEmbedding,TCell>::deserialize(
    Ref<TStream> rStream,
    Ref<DirichletPenaltyLHSIntegrand> rInstance,
    Ref<const TAllocator> rAllocator,
    tags::Binary) {
        using BS = cie::io::BinarySerializer;
        BS::deserialize(rStream, rInstance._ansatzSpace, rAllocator);
        BS::deserialize(rStream, rInstance._embedding, rAllocator);
        BS::deserialize(rStream, rInstance._cellInverseTransform, rAllocator);
        BS::deserialize(rStream, rInstance._cellJacobian, rAllocator);
        BS::deserialize(rStream, rInstance._penalty, rAllocator);
}


template <
    maths::Expression TDirichlet,
    maths::SliceableExpression TAnsatz,
    maths::Expression TEmbedding,
    CellLike TCell>
DirichletPenaltyRHSIntegrand<TDirichlet,TAnsatz,TEmbedding,TCell>::DirichletPenaltyRHSIntegrand() noexcept
    :   _ansatzSpace(),
        _embedding(),
        _cellInverseTransform(),
        _cellJacobian(),
        _dirichletFunctor(),
        _penalty(0)
{}


template <
    maths::Expression TDirichlet,
    maths::SliceableExpression TAnsatz,
    maths::Expression TEmbedding,
    CellLike TCell>
DirichletPenaltyRHSIntegrand<TDirichlet,TAnsatz,TEmbedding,TCell>::DirichletPenaltyRHSIntegrand(
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
    maths::SliceableExpression TAnsatz,
    maths::Expression TEmbedding,
    CellLike TCell>
void DirichletPenaltyRHSIntegrand<TDirichlet,TAnsatz,TEmbedding,TCell>::evaluate(
    ConstSpan in,
    Span out,
    BufferSpan buffer) const {
        Context context;
        const auto bounds = this->contextBounds();
        for (unsigned iDimension=0u; iDimension<ParametricDimension; ++iDimension) {
            context.data()[2 * iDimension] = 0;
            context.data()[2 * iDimension + 1] = bounds[iDimension];
        }
        const InContext inContext{context};
        return this->evaluate(in, out, buffer, inContext);
}


template <
    maths::Expression TDirichlet,
    maths::SliceableExpression TAnsatz,
    maths::Expression TEmbedding,
    CellLike TCell>
void DirichletPenaltyRHSIntegrand<TDirichlet,TAnsatz,TEmbedding,TCell>::evaluate(
    ConstSpan in,
    Span out,
    BufferSpan buffer,
    Ref<const InContext> rInContext) const {
        using AnsatzInContext = typename TAnsatz::InContext;
        const AnsatzInContext ansatzContext{rInContext.front()};
        const unsigned ansatzCount = _ansatzSpace.size(ansatzContext);
        const unsigned stateVariableCount = _dirichletFunctor.size();
        assert(out.size() == ansatzCount * stateVariableCount);
        assert(this->bufferSize(rInContext) <= buffer.size());

        Ptr<Value> pAnsatzBuffer = buffer.data();
        Ptr<Value> pEmbeddingBuffer = pAnsatzBuffer + ansatzCount;
        Ptr<Value> pCellParametricBuffer = pEmbeddingBuffer + TCell::PhysicalDimension;
        Ptr<Value> pDirichletBuffer = pCellParametricBuffer + TCell::ParametricDimension;
        Ptr<Value> pNestedBuffer = pDirichletBuffer + stateVariableCount;
        const Span ansatzBuffer(pAnsatzBuffer, ansatzCount);
        const Span boundaryTransformedBuffer(pEmbeddingBuffer, TCell::PhysicalDimension);
        const Span cellTransformedBuffer(pCellParametricBuffer, TCell::ParametricDimension);
        const Span dirichletBuffer(pDirichletBuffer, stateVariableCount);
        const BufferSpan nestedBuffer(pNestedBuffer, buffer.data() + buffer.size());

        const Value jacobian = std::abs(_cellJacobian.evaluateDeterminant(in, {}));
        const Value scaledPenalty = jacobian * _penalty;
        _embedding.evaluate(in, boundaryTransformedBuffer, nestedBuffer);
        _cellInverseTransform.evaluate(boundaryTransformedBuffer, cellTransformedBuffer, nestedBuffer);
        _ansatzSpace.evaluate(cellTransformedBuffer, ansatzBuffer, nestedBuffer, ansatzContext);
        _dirichletFunctor.evaluate(in, dirichletBuffer, nestedBuffer);

        for (unsigned iStateComponent=0u; iStateComponent<stateVariableCount; ++iStateComponent) {
            const Value dirichlet = dirichletBuffer[iStateComponent];
            std::transform(
                ansatzBuffer.begin(),
                ansatzBuffer.end(),
                out.begin() + iStateComponent * ansatzCount,
                [scaledPenalty, dirichlet] (const Value ansatz) -> Value {
                    return dirichlet * ansatz * scaledPenalty;
                });
        }
}


template <
    maths::Expression TDirichlet,
    maths::SliceableExpression TAnsatz,
    maths::Expression TEmbedding,
    CellLike TCell>
unsigned DirichletPenaltyRHSIntegrand<TDirichlet,TAnsatz,TEmbedding,TCell>::size() const noexcept {
    return _ansatzSpace.size() * _dirichletFunctor.size();
}


template <
    maths::Expression TDirichlet,
    maths::SliceableExpression TAnsatz,
    maths::Expression TEmbedding,
    CellLike TCell>
unsigned DirichletPenaltyRHSIntegrand<TDirichlet,TAnsatz,TEmbedding,TCell>::size(
    Ref<const InContext> rInContext) const noexcept {
        using AnsatzInContext = typename TAnsatz::InContext;
        return _ansatzSpace.size(AnsatzInContext{rInContext.front()}) * _dirichletFunctor.size();
}


template <
    maths::Expression TDirichlet,
    maths::SliceableExpression TAnsatz,
    maths::Expression TEmbedding,
    CellLike TCell>
unsigned DirichletPenaltyRHSIntegrand<TDirichlet,TAnsatz,TEmbedding,TCell>::bufferSize() const noexcept {
    const auto bounds = this->contextBounds();
    Context context;
    for (unsigned iDimension=0u; iDimension<ParametricDimension; ++iDimension) {
        context.data()[2 * iDimension] = 0;
        context.data()[2 * iDimension + 1] = bounds[iDimension];
    }
    const InContext inContext{context};
    return this->bufferSize(inContext);
}


template <
    maths::Expression TDirichlet,
    maths::SliceableExpression TAnsatz,
    maths::Expression TEmbedding,
    CellLike TCell>
unsigned DirichletPenaltyRHSIntegrand<TDirichlet,TAnsatz,TEmbedding,TCell>::bufferSize(
    Ref<const InContext> rInContext) const noexcept {
        using AnsatzInContext = typename TAnsatz::InContext;
        const AnsatzInContext ansatzContext{rInContext.front()};
        const unsigned ansatzCount = _ansatzSpace.size(ansatzContext);
        const unsigned stateVariableCount = _dirichletFunctor.size();
        const std::array<unsigned,4> nestedBufferRequirements {
            _dirichletFunctor.bufferSize(),
            _ansatzSpace.bufferSize(ansatzContext),
            _embedding.bufferSize(),
            _cellInverseTransform.bufferSize()};
        return ansatzCount
             + TCell::PhysicalDimension
             + TCell::ParametricDimension
             + stateVariableCount
             + *std::max_element(nestedBufferRequirements.begin(), nestedBufferRequirements.end());
}


template <
    maths::Expression TDirichlet,
    maths::SliceableExpression TAnsatz,
    maths::Expression TEmbedding,
    CellLike TCell>
std::array<
    typename DirichletPenaltyRHSIntegrand<TDirichlet,TAnsatz,TEmbedding,TCell>::Context::Value,
    DirichletPenaltyRHSIntegrand<TDirichlet,TAnsatz,TEmbedding,TCell>::ParametricDimension>
DirichletPenaltyRHSIntegrand<TDirichlet,TAnsatz,TEmbedding,TCell>::contextBounds() const noexcept {
    return _ansatzSpace.contextBounds();
}


template <
    maths::Expression TDirichlet,
    maths::SliceableExpression TAnsatz,
    maths::Expression TEmbedding,
    CellLike TCell>
typename DirichletPenaltyRHSIntegrand<TDirichlet,TAnsatz,TEmbedding,TCell>::OutContext
DirichletPenaltyRHSIntegrand<TDirichlet,TAnsatz,TEmbedding,TCell>::outContext(
    Ref<const InContext> rInContext) const noexcept {
        std::array<typename OutContext::Value,2 * (ParametricDimension + 1)> data{};
        for (unsigned iDimension=0u; iDimension<ParametricDimension; ++iDimension) {
            data[2 * iDimension] = rInContext.front().data()[2 * iDimension];
            data[2 * iDimension + 1] = rInContext.front().data()[2 * iDimension + 1];
        }
        data[2 * ParametricDimension] = 0;
        data[2 * ParametricDimension + 1] = _dirichletFunctor.size();
        return OutContext(data);
}


template <
    maths::Expression TDirichlet,
    maths::SliceableExpression TAnsatz,
    maths::Expression TEmbedding,
    CellLike TCell>
void DirichletPenaltyRHSIntegrand<TDirichlet,TAnsatz,TEmbedding,TCell>::serialize(
    Ref<cie::io::Traits::SerializerStream> rStream,
    tags::Binary) const {
        using BS = cie::io::BinarySerializer;
        BS::serialize(rStream, _ansatzSpace);
        BS::serialize(rStream, _embedding);
        BS::serialize(rStream, _cellInverseTransform);
        BS::serialize(rStream, _cellJacobian);
        BS::serialize(rStream, _dirichletFunctor);
        BS::serialize(rStream, _penalty);
}


template <
    maths::Expression TDirichlet,
    maths::SliceableExpression TAnsatz,
    maths::Expression TEmbedding,
    CellLike TCell>
template <class TStream, class TAllocator>
void DirichletPenaltyRHSIntegrand<TDirichlet,TAnsatz,TEmbedding,TCell>::deserialize(
    Ref<TStream> rStream,
    Ref<DirichletPenaltyRHSIntegrand> rInstance,
    Ref<const TAllocator> rAllocator,
    tags::Binary) {
        using BS = cie::io::BinarySerializer;
        BS::deserialize(rStream, rInstance._ansatzSpace, rAllocator);
        BS::deserialize(rStream, rInstance._embedding, rAllocator);
        BS::deserialize(rStream, rInstance._cellInverseTransform, rAllocator);
        BS::deserialize(rStream, rInstance._cellJacobian, rAllocator);
        BS::deserialize(rStream, rInstance._dirichletFunctor, rAllocator);
        BS::deserialize(rStream, rInstance._penalty, rAllocator);
}


} // namespace cie::fem
