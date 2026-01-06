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
{
    if constexpr (!maths::BufferedExpression<TAnsatzDerivatives>) {
        _ansatzDerivatives = nullptr;
    }
}


template <maths::Expression TAnsatzDerivatives>
LinearIsotropicStiffnessIntegrand<TAnsatzDerivatives>::LinearIsotropicStiffnessIntegrand(const Value modulus,
                                                                                         RightRef<TAnsatzDerivatives> rAnsatzDerivatives)
requires maths::BufferedExpression<TAnsatzDerivatives>
    : _modulus(modulus),
      _ansatzDerivatives(std::move(rAnsatzDerivatives)),
      _buffer()
{
}


template <maths::Expression TAnsatzDerivatives>
LinearIsotropicStiffnessIntegrand<TAnsatzDerivatives>::LinearIsotropicStiffnessIntegrand(const Value modulus,
                                                                                         RightRef<TAnsatzDerivatives> rAnsatzDerivatives,
                                                                                         std::span<Value> buffer)
requires maths::BufferedExpression<TAnsatzDerivatives>
    : LinearIsotropicStiffnessIntegrand(modulus, std::move(rAnsatzDerivatives))
{
    this->setBuffer(buffer);
}


template <maths::Expression TAnsatzDerivatives>
LinearIsotropicStiffnessIntegrand<TAnsatzDerivatives>::LinearIsotropicStiffnessIntegrand(const Value modulus,
                                                                                         Ref<const TAnsatzDerivatives> rAnsatzDerivatives)
requires (!maths::BufferedExpression<TAnsatzDerivatives>)
    : _modulus(modulus),
      _ansatzDerivatives(&rAnsatzDerivatives),
      _buffer()
{
}


template <maths::Expression TAnsatzDerivatives>
LinearIsotropicStiffnessIntegrand<TAnsatzDerivatives>::LinearIsotropicStiffnessIntegrand(const Value modulus,
                                                                                         Ref<const TAnsatzDerivatives> rAnsatzDerivatives,
                                                                                         std::span<Value> buffer)
requires (!maths::BufferedExpression<TAnsatzDerivatives>)
    : LinearIsotropicStiffnessIntegrand(modulus, rAnsatzDerivatives)
{
    this->setBuffer(buffer);
}


template <maths::Expression TAnsatzDerivatives>
void LinearIsotropicStiffnessIntegrand<TAnsatzDerivatives>::evaluate(ConstSpan in, Span out) const {
    CIE_OUT_OF_RANGE_CHECK(this->getMinBufferSize() <= _buffer.size())

    if constexpr (!maths::BufferedExpression<TAnsatzDerivatives>) {
        CIE_CHECK_POINTER(_ansatzDerivatives)
    }

    Ref<const TAnsatzDerivatives> rAnsatzDerivatives = utils::getRef(_ansatzDerivatives);
    const unsigned derivativeComponentCount = rAnsatzDerivatives.size();
    const unsigned ansatzCount = derivativeComponentCount / Dimension;
    rAnsatzDerivatives.evaluate(in, {_buffer.data(), _buffer.data() + rAnsatzDerivatives.size()});

    using EigenDenseMatrix = Eigen::Matrix<Value,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>;
    using EigenAdaptor = Eigen::Map<EigenDenseMatrix>;

    EigenAdaptor derivativeAdaptor(_buffer.data(), Dimension, ansatzCount);
    EigenAdaptor outputAdaptor(out.data(), ansatzCount, ansatzCount);

    outputAdaptor = derivativeAdaptor.transpose() * _modulus * derivativeAdaptor;
}


template <maths::Expression TAnsatzDerivatives>
unsigned LinearIsotropicStiffnessIntegrand<TAnsatzDerivatives>::size() const {
    const auto derivativeComponentCount = utils::getRef(_ansatzDerivatives).size();
    const auto ansatzCount = derivativeComponentCount / Dimension;
    return ansatzCount * ansatzCount;
}


template <maths::Expression TAnsatzDerivatives>
void LinearIsotropicStiffnessIntegrand<TAnsatzDerivatives>::setBuffer(std::span<Value> buffer) {
    CIE_OUT_OF_RANGE_CHECK(this->getMinBufferSize() <= buffer.size())
    if constexpr (maths::BufferedExpression<TAnsatzDerivatives>) {
        unsigned ansatzBufferSize = utils::getRef(_ansatzDerivatives).getMinBufferSize();
        utils::getRef(_ansatzDerivatives).setBuffer({
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
unsigned LinearIsotropicStiffnessIntegrand<TAnsatzDerivatives>::getMinBufferSize() const noexcept {
    unsigned bufferSize = utils::getRef(_ansatzDerivatives).size();
    if constexpr (maths::BufferedExpression<TAnsatzDerivatives>) {
        bufferSize += utils::getRef(_ansatzDerivatives).getMinBufferSize();
    }
    return bufferSize;
}


} // namespace cie::fem
