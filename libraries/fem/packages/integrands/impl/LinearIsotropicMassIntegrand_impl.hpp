#pragma once

// --- External Includes ---
#include <Eigen/Dense> // Eigen::Map

// help the language server
#include "packages/integrands/inc/LinearIsotropicMassIntegrand.hpp"


namespace cie::fem {


template <maths::Expression TAnsatzSpace>
LinearIsotropicMassIntegrand<TAnsatzSpace>::LinearIsotropicMassIntegrand()
    : LinearIsotropicMassIntegrand(0, nullptr)
{
}


template <maths::Expression TAnsatzSpace>
LinearIsotropicMassIntegrand<TAnsatzSpace>::LinearIsotropicMassIntegrand(const Value modulus,
                                                                         Ref<const TAnsatzSpace> rAnsatzSpace)
    : _modulus(modulus),
      _pAnsatzSpace(&rAnsatzSpace),
      _buffer()
{
}


template <maths::Expression TAnsatzSpace>
LinearIsotropicMassIntegrand<TAnsatzSpace>::LinearIsotropicMassIntegrand(const Value modulus,
                                                                         Ref<const TAnsatzSpace> rAnsatzSpace,
                                                                         std::span<Value> buffer)
    : LinearIsotropicMassIntegrand(modulus, rAnsatzSpace)
{
    this->setBuffer(buffer);
}


template <maths::Expression TAnsatzSpace>
void LinearIsotropicMassIntegrand<TAnsatzSpace>::evaluate(ConstSpan in, Span out) const
{
    CIE_OUT_OF_RANGE_CHECK(this->getMinBufferSize() <= _buffer.size())
    CIE_CHECK_POINTER(_pAnsatzSpace)

    Ref<const TAnsatzSpace> rAnsatzSpace = *_pAnsatzSpace;
    const unsigned ansatzCount = rAnsatzSpace.size();
    rAnsatzSpace.evaluate(in, {_buffer.data(), _buffer.data() + _pAnsatzSpace->size()});

    using EigenDenseMatrix = Eigen::Matrix<Value,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>;
    using EigenAdaptor = Eigen::Map<EigenDenseMatrix>;

    EigenAdaptor ansatzAdaptor(_buffer.data(), ansatzCount, 1);
    EigenAdaptor outputAdaptor(out.data(), ansatzCount, ansatzCount);

    outputAdaptor = ansatzAdaptor * _modulus * ansatzAdaptor.transpose();
}


template <maths::Expression TAnsatzSpace>
unsigned LinearIsotropicMassIntegrand<TAnsatzSpace>::size() const
{
    const auto ansatzCount = _pAnsatzSpace->size();
    return ansatzCount * ansatzCount;
}


template <maths::Expression TAnsatzSpace>
void LinearIsotropicMassIntegrand<TAnsatzSpace>::setBuffer(std::span<Value> buffer)
{
    CIE_OUT_OF_RANGE_CHECK(this->getMinBufferSize() <= buffer.size())
    _buffer = buffer;
}


template <maths::Expression TAnsatzSpace>
unsigned LinearIsotropicMassIntegrand<TAnsatzSpace>::getMinBufferSize() const noexcept
{
    return _pAnsatzSpace->size();
}


} // namespace cie::fem
