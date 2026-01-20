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
      _ansatzDerivatives(),
      _buffer()
{}


template <maths::Expression TAnsatzDerivatives>
LinearIsotropicStiffnessIntegrand<TAnsatzDerivatives>::LinearIsotropicStiffnessIntegrand(
    const Value modulus,
    RightRef<TAnsatzDerivatives> rAnsatzDerivatives) noexcept
    : _modulus(modulus),
      _ansatzDerivatives(std::move(rAnsatzDerivatives)),
      _buffer()
{}


template <maths::Expression TAnsatzDerivatives>
LinearIsotropicStiffnessIntegrand<TAnsatzDerivatives>::LinearIsotropicStiffnessIntegrand(
    const Value modulus,
    RightRef<TAnsatzDerivatives> rAnsatzDerivatives,
    std::span<Value> buffer)
requires (isBuffered)
    : LinearIsotropicStiffnessIntegrand(modulus, std::move(rAnsatzDerivatives))
{
    this->setBuffer(buffer);
}


template <maths::Expression TAnsatzDerivatives>
LinearIsotropicStiffnessIntegrand<TAnsatzDerivatives>::LinearIsotropicStiffnessIntegrand(const Value modulus,
                                                                                         Ref<const TAnsatzDerivatives> rAnsatzDerivatives)
    : _modulus(modulus),
      _ansatzDerivatives(rAnsatzDerivatives),
      _buffer()
{
}


template <maths::Expression TAnsatzDerivatives>
LinearIsotropicStiffnessIntegrand<TAnsatzDerivatives>::LinearIsotropicStiffnessIntegrand(
    const Value modulus,
    Ref<const TAnsatzDerivatives> rAnsatzDerivatives,
    std::span<Value> buffer)
requires (isBuffered)
    : LinearIsotropicStiffnessIntegrand(modulus, rAnsatzDerivatives)
{
    this->setBuffer(buffer);
}


template <maths::Expression TAnsatzDerivatives>
void LinearIsotropicStiffnessIntegrand<TAnsatzDerivatives>::evaluate(ConstSpan in, Span out) const {

    unsigned derivativeComponentCount = 0u;
    if constexpr (maths::StaticExpression<TAnsatzDerivatives>) {
        derivativeComponentCount = TAnsatzDerivatives::size();
    } else {
        derivativeComponentCount = _ansatzDerivatives.size();
    }

    const unsigned ansatzCount = derivativeComponentCount / Dimension;
    _ansatzDerivatives.evaluate(in, {_buffer.data(), derivativeComponentCount});

    using EigenDenseMatrix = Eigen::Matrix<Value,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>;
    using EigenAdaptor = Eigen::Map<EigenDenseMatrix>;

    EigenAdaptor derivativeAdaptor(
        _buffer.data(),
        Dimension,
        ansatzCount);
    EigenAdaptor outputAdaptor(
        out.data(),
        ansatzCount,
        ansatzCount);

    //derivativeTransposeAdaptor.transposeInPlace();
    outputAdaptor.noalias() = _modulus * derivativeAdaptor.transpose().lazyProduct(derivativeAdaptor);
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
void LinearIsotropicStiffnessIntegrand<TAnsatzDerivatives>::setBuffer(std::span<Value> buffer)
requires (isBuffered) {
    CIE_OUT_OF_RANGE_CHECK(this->getMinBufferSize() <= buffer.size())
    if constexpr (maths::BufferedExpression<TAnsatzDerivatives>) {
        unsigned ansatzBufferSize = _ansatzDerivatives.getMinBufferSize();
        _ansatzDerivatives.setBuffer({
            buffer.data(),
            ansatzBufferSize});
        _buffer = {
            buffer.data() + ansatzBufferSize,
            buffer.data() + buffer.size()};
    } else {
        _buffer = buffer;
    }
}


template <maths::Expression TAnsatzDerivatives>
unsigned LinearIsotropicStiffnessIntegrand<TAnsatzDerivatives>::getMinBufferSize() const noexcept
requires (isBuffered) {
    unsigned bufferSize = _ansatzDerivatives.size();
    if constexpr (maths::BufferedExpression<TAnsatzDerivatives>) {
        bufferSize += _ansatzDerivatives.getMinBufferSize();
    }
    return bufferSize;
}


} // namespace cie::fem
