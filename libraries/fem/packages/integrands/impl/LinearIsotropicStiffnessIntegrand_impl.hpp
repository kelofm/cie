#pragma once

// --- External Includes ---
#include <Eigen/Dense> // Eigen::Map

// --- FEM Includes ---
#include "packages/integrands/inc/LinearIsotropicStiffnessIntegrand.hpp"

// --- Utility Includes ---
#include "packages/io/inc/Serializer.hpp"


namespace cie::fem {


template <maths::Expression TAD, maths::SpatialTransform TT>
LinearIsotropicStiffnessIntegrand<TAD,TT>::LinearIsotropicStiffnessIntegrand()
    :   _ansatzDerivatives(),
        _jacobian(),
        _jacobianInverse(),
        _modulus(0)
{}


template <maths::Expression TAD, maths::SpatialTransform TT>
LinearIsotropicStiffnessIntegrand<TAD,TT>::LinearIsotropicStiffnessIntegrand(
    const Value modulus,
    RightRef<TAD> rAnsatzDerivatives,
    RightRef<Jacobian> rJacobian,
    RightRef<JacobianInverse> rJacobianInverse) noexcept
    :   _ansatzDerivatives(std::move(rAnsatzDerivatives)),
        _jacobian(std::move(rJacobian)),
        _jacobianInverse(std::move(rJacobianInverse)),
        _modulus(modulus)
{}


template <maths::Expression TAD, maths::SpatialTransform TT>
LinearIsotropicStiffnessIntegrand<TAD,TT>::LinearIsotropicStiffnessIntegrand(
    const Value modulus,
    Ref<const TAD> rAnsatzDerivatives,
    Ref<const Jacobian> rJacobian,
    Ref<const JacobianInverse> rJacobianInverse)
        : LinearIsotropicStiffnessIntegrand(
            modulus,
            TAD(rAnsatzDerivatives),
            Jacobian(rJacobian),
            JacobianInverse(rJacobianInverse))
{}


/// @name Common Expression Interface
/// @{


template <maths::Expression TAD, maths::SpatialTransform TT>
void LinearIsotropicStiffnessIntegrand<TAD,TT>::evaluate(
    ConstSpan in,
    Span out,
    BufferSpan buffer) const {
        const typename Context::Value ansatzSetSize = _ansatzDerivatives.ansatzSet().size();
        Context context;
        for (auto it=context.data().begin(); it<context.data().end(); it+=2) {
            *it = 0;
            *(it + 1) = ansatzSetSize;
        }
        const InContext inContext{context,context};
        return this->evaluate(in, out, buffer, inContext);
}


template <maths::Expression TAD, maths::SpatialTransform TT>
void LinearIsotropicStiffnessIntegrand<TAD,TT>::evaluate(
    ConstSpan in,
    Span out,
    BufferSpan buffer,
    Ref<const Context> rContext) const {
        const InContext inContext{rContext,rContext};
        return this->evaluate(in, out, buffer, inContext);
}


template <maths::Expression TAD, maths::SpatialTransform TT>
void LinearIsotropicStiffnessIntegrand<TAD,TT>::evaluate(
    ConstSpan in,
    Span out,
    BufferSpan buffer,
    Ref<const InContext> rInContext) const {
        using DerivativeInContext = typename TAD::InContext;
        const DerivativeInContext rowContext{rInContext[0]};
        const DerivativeInContext columnContext{rInContext[1]};
        const std::size_t rowDerivativeCount = _ansatzDerivatives.size(rowContext);
        const std::size_t columnDerivativeCount = _ansatzDerivatives.size(columnContext);
        assert(out.size() == this->size(rInContext));
        assert(this->bufferSize(rInContext) <= buffer.size());
        Span rowDerivativeBuffer(buffer.data(), rowDerivativeCount);
        Span columnDerivativeBuffer(buffer.data() + rowDerivativeCount, columnDerivativeCount);
        Span nestedBuffer(
            buffer.data() + rowDerivativeCount + columnDerivativeCount,
            buffer.data() + buffer.size());

        // Evaluate the jacobian's inverse.
        // \frac{d\xi_i}{dx_j}
        static_assert(JacobianInverse::size() == ParametricDimension * ParametricDimension);
        std::array<Value,JacobianInverse::size()> jacobianInverse;
        _jacobianInverse.evaluate(
            in,
            jacobianInverse,
            {});
        Eigen::Map<Eigen::Matrix<Value,ParametricDimension,ParametricDimension,Eigen::RowMajor>> jacobianInverseAdaptor(
            jacobianInverse.data());

        // Compute \frac{d\xi_j}{dx_i} \frac{d\xi_j}{dx_k}
        {
            std::array<Value,JacobianInverse::size()> productBuffer;
            Eigen::Map<Eigen::Matrix<Value,ParametricDimension,ParametricDimension,Eigen::RowMajor>> productBufferAdaptor(
                productBuffer.data());
            productBufferAdaptor.noalias() = jacobianInverseAdaptor.transpose().lazyProduct(jacobianInverseAdaptor);
            jacobianInverse = productBuffer;
        }

        // Compute the jacobian's determinant.
        // det(\frac{dx_i}{d\xi_j})
        const Value jacobianDeterminant = std::abs(_jacobian.evaluateDeterminant(
            in,
            {}));

        const auto rowAnsatzCount = rowDerivativeCount / ParametricDimension;
        const auto columnAnsatzCount = columnDerivativeCount / ParametricDimension;

        // Evaluate the ansatz functions' derivatives.
        // \frac{dN_i}{d\xi_j}
        _ansatzDerivatives.evaluate(
            in,
            rowDerivativeBuffer,
            nestedBuffer,
            rowContext);
        _ansatzDerivatives.evaluate(
            in,
            columnDerivativeBuffer,
            nestedBuffer,
            columnContext);

        using EigenDenseMatrix = Eigen::Matrix<Value,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>;
        using EigenAdaptor = Eigen::Map<EigenDenseMatrix>;
        EigenAdaptor derivativeAdaptor(
            rowDerivativeBuffer.data(),
            ParametricDimension,
            rowAnsatzCount);
        EigenAdaptor columnDerivativeAdaptor(
            columnDerivativeBuffer.data(),
            ParametricDimension,
            columnAnsatzCount);

        // Compute the integrand.
        // det(\frac{dx_i}{d\xi_j}) m \frac{dN_j}{dx_i} \frac{dx_j}{d\xi_k} \frac{dN_k}{d\xi_l}
        EigenAdaptor outputAdaptor(
            out.data(),
            rowAnsatzCount,
            columnAnsatzCount);
        outputAdaptor.noalias() = (jacobianDeterminant * _modulus) * derivativeAdaptor
            .transpose()
            .lazyProduct(
                jacobianInverseAdaptor.lazyProduct(
                    columnDerivativeAdaptor));
}


/// @}
/// @name Runtime Expression Interface
/// @{


template <maths::Expression TAD, maths::SpatialTransform TT>
unsigned LinearIsotropicStiffnessIntegrand<TAD,TT>::size() const noexcept
requires (!maths::StaticExpression<TAD>) {
    const auto derivativeComponentCount = _ansatzDerivatives.size();
    const auto ansatzCount = derivativeComponentCount / ParametricDimension;
    return ansatzCount * ansatzCount;
}


template <maths::Expression TAD, maths::SpatialTransform TT>
unsigned LinearIsotropicStiffnessIntegrand<TAD,TT>::bufferSize() const noexcept
requires (!maths::StaticExpression<TAD>) {
    return 2 * _ansatzDerivatives.size() + _ansatzDerivatives.bufferSize();
}


template <maths::Expression TAD, maths::SpatialTransform TT>
unsigned LinearIsotropicStiffnessIntegrand<TAD,TT>::size(Ref<const Context> rContext) const noexcept
requires (!maths::StaticExpression<TAD>) {
    const InContext inContext{rContext,rContext};
    return this->size(inContext);
}


template <maths::Expression TAD, maths::SpatialTransform TT>
unsigned LinearIsotropicStiffnessIntegrand<TAD,TT>::size(Ref<const InContext> rInContext) const noexcept
requires (!maths::StaticExpression<TAD>) {
    using DerivativeInContext = typename TAD::InContext;
    const unsigned rowCount = _ansatzDerivatives.size(DerivativeInContext{rInContext[0]}) / ParametricDimension;
    const unsigned columnCount = _ansatzDerivatives.size(DerivativeInContext{rInContext[1]}) / ParametricDimension;
    return rowCount * columnCount;
}


template <maths::Expression TAD, maths::SpatialTransform TT>
unsigned LinearIsotropicStiffnessIntegrand<TAD,TT>::bufferSize(Ref<const Context> rContext) const noexcept
requires (!maths::StaticExpression<TAD>) {
    const InContext inContext{rContext,rContext};
    return this->bufferSize(inContext);
}


template <maths::Expression TAD, maths::SpatialTransform TT>
unsigned LinearIsotropicStiffnessIntegrand<TAD,TT>::bufferSize(Ref<const InContext> rInContext) const noexcept
requires (!maths::StaticExpression<TAD>) {
    using DerivativeInContext = typename TAD::InContext;
    const DerivativeInContext rowContext{rInContext[0]};
    const DerivativeInContext columnContext{rInContext[1]};
    const unsigned rowSize = _ansatzDerivatives.size(rowContext);
    const unsigned columnSize = _ansatzDerivatives.size(columnContext);
    return rowSize + columnSize + std::max(
        _ansatzDerivatives.bufferSize(rowContext),
        _ansatzDerivatives.bufferSize(columnContext));
}


template <maths::Expression TAD, maths::SpatialTransform TT>
std::array<
    typename LinearIsotropicStiffnessIntegrand<TAD,TT>::Context::Value,
    LinearIsotropicStiffnessIntegrand<TAD,TT>::ParametricDimension>
LinearIsotropicStiffnessIntegrand<TAD,TT>::contextBounds() const noexcept
requires (!maths::StaticExpression<TAD>) {
    return _ansatzDerivatives.contextBounds();
}


template <maths::Expression TAD, maths::SpatialTransform TT>
typename LinearIsotropicStiffnessIntegrand<TAD,TT>::OutContext
LinearIsotropicStiffnessIntegrand<TAD,TT>::outContext(Ref<const InContext> rInContext) const noexcept
requires (!maths::StaticExpression<TAD>) {
    std::array<typename OutContext::Value,4 * ParametricDimension> data{};
    for (unsigned iDimension=0u; iDimension<ParametricDimension; ++iDimension) {
        // Output is row-major: column indices are the fastest-varying coordinates.
        data[2 * iDimension] = rInContext[1].data()[2 * iDimension];
        data[2 * iDimension + 1] = rInContext[1].data()[2 * iDimension + 1];
        data[2 * (ParametricDimension + iDimension)] = rInContext[0].data()[2 * iDimension];
        data[2 * (ParametricDimension + iDimension) + 1] = rInContext[0].data()[2 * iDimension + 1];
    }
    return OutContext(data);
}


/// @}
/// @name Static Expression Interface
/// @{


template <maths::Expression TAD, maths::SpatialTransform TT>
constexpr unsigned LinearIsotropicStiffnessIntegrand<TAD,TT>::size() noexcept
requires (maths::StaticExpression<TAD>) {
    const auto derivativeComponentCount = TAD::size();
    const auto ansatzCount = derivativeComponentCount / ParametricDimension;
    return ansatzCount * ansatzCount;
}


template <maths::Expression TAD, maths::SpatialTransform TT>
constexpr unsigned LinearIsotropicStiffnessIntegrand<TAD,TT>::bufferSize() noexcept
requires (maths::StaticExpression<TAD>) {
    return 2 * TAD::size() + TAD::bufferSize();
}


template <maths::Expression TAD, maths::SpatialTransform TT>
constexpr unsigned LinearIsotropicStiffnessIntegrand<TAD,TT>::size(Ref<const Context> rContext) noexcept
requires (maths::StaticExpression<TAD>) {
    const InContext inContext{rContext,rContext};
    return size(inContext);
}


template <maths::Expression TAD, maths::SpatialTransform TT>
constexpr unsigned LinearIsotropicStiffnessIntegrand<TAD,TT>::size(Ref<const InContext> rInContext) noexcept
requires (maths::StaticExpression<TAD>) {
    using DerivativeInContext = typename TAD::InContext;
    const unsigned rowCount = TAD::size(DerivativeInContext{rInContext[0]}) / ParametricDimension;
    const unsigned columnCount = TAD::size(DerivativeInContext{rInContext[1]}) / ParametricDimension;
    return rowCount * columnCount;
}


template <maths::Expression TAD, maths::SpatialTransform TT>
constexpr unsigned LinearIsotropicStiffnessIntegrand<TAD,TT>::bufferSize(Ref<const Context> rContext) noexcept
requires (maths::StaticExpression<TAD>) {
    const InContext inContext{rContext,rContext};
    return bufferSize(inContext);
}


template <maths::Expression TAD, maths::SpatialTransform TT>
constexpr unsigned LinearIsotropicStiffnessIntegrand<TAD,TT>::bufferSize(Ref<const InContext> rInContext) noexcept
requires (maths::StaticExpression<TAD>) {
    using DerivativeInContext = typename TAD::InContext;
    const DerivativeInContext rowContext{rInContext[0]};
    const DerivativeInContext columnContext{rInContext[1]};
    const unsigned rowSize = TAD::size(rowContext);
    const unsigned columnSize = TAD::size(columnContext);
    return rowSize + columnSize + std::max(
        TAD::bufferSize(rowContext),
        TAD::bufferSize(columnContext));
}


template <maths::Expression TAD, maths::SpatialTransform TT>
constexpr std::array<
    typename LinearIsotropicStiffnessIntegrand<TAD,TT>::Context::Value,
    LinearIsotropicStiffnessIntegrand<TAD,TT>::ParametricDimension>
LinearIsotropicStiffnessIntegrand<TAD,TT>::contextBounds() noexcept
requires (maths::StaticExpression<TAD>) {
    return TAD::contextBounds();
}


template <maths::Expression TAD, maths::SpatialTransform TT>
constexpr typename LinearIsotropicStiffnessIntegrand<TAD,TT>::OutContext
LinearIsotropicStiffnessIntegrand<TAD,TT>::outContext(Ref<const InContext> rInContext) noexcept
requires (maths::StaticExpression<TAD>) {
    std::array<typename OutContext::Value,4 * ParametricDimension> data{};
    for (unsigned iDimension=0u; iDimension<ParametricDimension; ++iDimension) {
        data[2 * iDimension] = rInContext[1].data()[2 * iDimension];
        data[2 * iDimension + 1] = rInContext[1].data()[2 * iDimension + 1];
        data[2 * (ParametricDimension + iDimension)] = rInContext[0].data()[2 * iDimension];
        data[2 * (ParametricDimension + iDimension) + 1] = rInContext[0].data()[2 * iDimension + 1];
    }
    return OutContext(data);
}


/// @}
/// @name Serialization
/// @{


template <maths::Expression TAD, maths::SpatialTransform TT>
void LinearIsotropicStiffnessIntegrand<TAD,TT>::serialize(
    Ref<cie::io::Traits::SerializerStream> rStream,
    tags::Binary) const {
        using BS = cie::io::BinarySerializer;
        BS::serialize(
            rStream,
            _ansatzDerivatives);
        BS::serialize(
            rStream,
            _jacobian);
        BS::serialize(
            rStream,
            _jacobianInverse);
        BS::serialize(
            rStream,
            _modulus);
}


template <maths::Expression TAD, maths::SpatialTransform TT>
template <class TStream, class TAllocator>
void LinearIsotropicStiffnessIntegrand<TAD,TT>::deserialize(
    Ref<TStream> rStream,
    Ref<LinearIsotropicStiffnessIntegrand> rInstance,
    Ref<const TAllocator> rAllocator,
    tags::Binary) {
        using BS = cie::io::BinarySerializer;
        BS::deserialize(
            rStream,
            rInstance._ansatzDerivatives,
            rAllocator);
        BS::deserialize(
            rStream,
            rInstance._jacobian,
            rAllocator);
        BS::deserialize(
            rStream,
            rInstance._jacobianInverse,
            rAllocator);
        BS::deserialize(
            rStream,
            rInstance._modulus,
            rAllocator);
}


/// @}


} // namespace cie::fem
