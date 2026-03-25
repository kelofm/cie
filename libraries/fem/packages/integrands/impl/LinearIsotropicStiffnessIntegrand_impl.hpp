#pragma once

// --- External Includes ---
#include <Eigen/Dense> // Eigen::Map

// help the language server
#include "packages/integrands/inc/LinearIsotropicStiffnessIntegrand.hpp"

// --- Utility Includes ---
#include "packages/stl_extension/inc/getReference.hpp"


namespace cie::fem {


template <maths::Expression TAnsatzDerivatives>
LinearIsotropicStiffnessIntegrand<TAnsatzDerivatives>::LinearIsotropicStiffnessIntegrand()
    : _modulus(0),
      _ansatzDerivatives()
{}


template <maths::Expression TAnsatzDerivatives>
LinearIsotropicStiffnessIntegrand<TAnsatzDerivatives>::LinearIsotropicStiffnessIntegrand(
    const Value modulus,
    RightRef<TAnsatzDerivatives> rAnsatzDerivatives) noexcept
    : _modulus(modulus),
      _ansatzDerivatives(std::move(rAnsatzDerivatives))
{}


template <maths::Expression TAnsatzDerivatives>
LinearIsotropicStiffnessIntegrand<TAnsatzDerivatives>::LinearIsotropicStiffnessIntegrand(const Value modulus,
                                                                                         Ref<const TAnsatzDerivatives> rAnsatzDerivatives)
    : _modulus(modulus),
      _ansatzDerivatives(rAnsatzDerivatives)
{}


template <maths::Expression TAnsatzDerivatives>
void LinearIsotropicStiffnessIntegrand<TAnsatzDerivatives>::evaluate(
    ConstSpan in,
    Span out,
    BufferSpan buffer) const {
        Span ansatzDerivativeBuffer, nestedBuffer;
        [[maybe_unused]] std::array<Value,maths::StaticExpressionSize<TAnsatzDerivatives>::size> ansatzDerivativeBufferArray;
        [[maybe_unused]] std::array<Value,maths::StaticExpressionSize<TAnsatzDerivatives>::bufferSize> nestedBufferArray;

        if constexpr (maths::StaticExpression<TAnsatzDerivatives>) {
            ansatzDerivativeBuffer = ansatzDerivativeBufferArray;
            nestedBuffer = nestedBufferArray;
        } else {
            const unsigned derivativeComponentCount = _ansatzDerivatives.size();
            ansatzDerivativeBuffer = Span(
                buffer.data(),
                derivativeComponentCount);
            nestedBuffer = Span(
                buffer.data() + derivativeComponentCount,
                buffer.data() + buffer.size());
        }

        if constexpr (maths::StaticExpression<TAnsatzDerivatives>) {
            constexpr unsigned derivativeComponentCount = TAnsatzDerivatives::size();
            constexpr unsigned ansatzCount = derivativeComponentCount / Dimension;
            _ansatzDerivatives.evaluate(in, ansatzDerivativeBuffer, nestedBuffer);

            using DerivativeMatrix = Eigen::Matrix<Value,Dimension,ansatzCount,Eigen::RowMajor>;
            using OutputMatrix = Eigen::Matrix<Value,ansatzCount,ansatzCount,Eigen::RowMajor>;
            Eigen::Map<DerivativeMatrix> derivativeAdaptor(ansatzDerivativeBuffer.data());
            Eigen::Map<OutputMatrix> outputAdaptor(
                out.data(),
                ansatzCount,
                ansatzCount);
            outputAdaptor.noalias() = _modulus * derivativeAdaptor.transpose().lazyProduct(derivativeAdaptor);
        } else {
            const unsigned derivativeComponentCount = ansatzDerivativeBuffer.size();
            const unsigned ansatzCount = derivativeComponentCount / Dimension;
            _ansatzDerivatives.evaluate(in, ansatzDerivativeBuffer, nestedBuffer);

            using EigenDenseMatrix = Eigen::Matrix<Value,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>;
            using EigenAdaptor = Eigen::Map<EigenDenseMatrix>;
            EigenAdaptor derivativeAdaptor(
                ansatzDerivativeBuffer.data(),
                Dimension,
                ansatzCount);
            EigenAdaptor outputAdaptor(
                out.data(),
                ansatzCount,
                ansatzCount);
            outputAdaptor.noalias() = _modulus * derivativeAdaptor.transpose().lazyProduct(derivativeAdaptor);
        }
}


template <maths::Expression TAnsatzDerivatives>
unsigned LinearIsotropicStiffnessIntegrand<TAnsatzDerivatives>::size() const noexcept
requires (!maths::StaticExpression<TAnsatzDerivatives>) {
    const auto derivativeComponentCount = _ansatzDerivatives.size();
    const auto ansatzCount = derivativeComponentCount / Dimension;
    return ansatzCount * ansatzCount;
}


template <maths::Expression TAnsatzDerivatives>
constexpr unsigned LinearIsotropicStiffnessIntegrand<TAnsatzDerivatives>::size() noexcept
requires (maths::StaticExpression<TAnsatzDerivatives>) {
    const auto derivativeComponentCount = TAnsatzDerivatives::size();
    const auto ansatzCount = derivativeComponentCount / Dimension;
    return ansatzCount * ansatzCount;
}


template <maths::Expression TAnsatzDerivatives>
constexpr unsigned LinearIsotropicStiffnessIntegrand<TAnsatzDerivatives>::bufferSize() noexcept
requires (maths::StaticExpression<TAnsatzDerivatives>) {
    return TAnsatzDerivatives::bufferSize();
}


template <maths::Expression TAnsatzDerivatives>
unsigned LinearIsotropicStiffnessIntegrand<TAnsatzDerivatives>::bufferSize() const noexcept
requires (!maths::StaticExpression<TAnsatzDerivatives>) {
    return _ansatzDerivatives.size() + _ansatzDerivatives.bufferSize();
}


} // namespace cie::fem
