#pragma once

// --- External Includes ---
#include <Eigen/Dense> // Eigen::Map

// help the language server
#include "packages/integrands/inc/LinearIsotropicStiffnessIntegrand.hpp"


namespace cie::fem {


template <maths::Expression TAD, maths::SpatialTransform TT>
LinearIsotropicStiffnessIntegrand<TAD,TT>::LinearIsotropicStiffnessIntegrand()
    :   _modulus(0),
        _ansatzDerivatives(),
        _jacobian(),
        _jacobianInverse()
{}


template <maths::Expression TAD, maths::SpatialTransform TT>
LinearIsotropicStiffnessIntegrand<TAD,TT>::LinearIsotropicStiffnessIntegrand(
    const Value modulus,
    RightRef<TAD> rAnsatzDerivatives,
    RightRef<Jacobian> rJacobian,
    RightRef<JacobianInverse> rJacobianInverse) noexcept
    :   _modulus(modulus),
        _ansatzDerivatives(std::move(rAnsatzDerivatives)),
        _jacobian(std::move(rJacobian)),
        _jacobianInverse(std::move(rJacobianInverse))
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


template <maths::Expression TAD, maths::SpatialTransform TT>
void LinearIsotropicStiffnessIntegrand<TAD,TT>::evaluate(
    ConstSpan in,
    Span out,
    BufferSpan buffer) const {
        Span ansatzDerivativeBuffer, nestedBuffer;
        [[maybe_unused]] std::array<Value,maths::StaticExpressionSize<TAD>::size> ansatzDerivativeBufferArray;
        [[maybe_unused]] std::array<Value,maths::StaticExpressionSize<TAD>::bufferSize> nestedBufferArray;

        if constexpr (maths::StaticExpression<TAD>) {
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

        std::array<Value,JacobianInverse::size()> jacobianInverse;
        _jacobianInverse.evaluate(
            in,
            jacobianInverse,
            {});
        Eigen::Map<Eigen::Matrix<Value,Dimension,Dimension,Eigen::RowMajor>> jacobianInverseAdaptor(
            jacobianInverse.data());

        {
            std::array<Value,JacobianInverse::size()> productBuffer;
            Eigen::Map<Eigen::Matrix<Value,Dimension,Dimension,Eigen::RowMajor>> productBufferAdaptor(
                productBuffer.data());
            productBufferAdaptor.noalias() = jacobianInverseAdaptor.transpose().lazyProduct(jacobianInverseAdaptor);
            jacobianInverse = productBuffer;
        }

        const Value jacobianDeterminant = std::abs(_jacobian.evaluateDeterminant(
            in,
            {}));

        if constexpr (maths::StaticExpression<TAD>) {
            constexpr unsigned derivativeComponentCount = TAD::size();
            constexpr unsigned ansatzCount = derivativeComponentCount / Dimension;
            _ansatzDerivatives.evaluate(in, ansatzDerivativeBuffer, nestedBuffer);

            using DerivativeMatrix = Eigen::Matrix<Value,Dimension,ansatzCount,Eigen::RowMajor>;
            using OutputMatrix = Eigen::Matrix<Value,ansatzCount,ansatzCount,Eigen::RowMajor>;
            Eigen::Map<DerivativeMatrix> derivativeAdaptor(ansatzDerivativeBuffer.data());
            Eigen::Map<OutputMatrix> outputAdaptor(
                out.data(),
                ansatzCount,
                ansatzCount);
            outputAdaptor.noalias() = (jacobianDeterminant * _modulus) * derivativeAdaptor
                .transpose()
                .lazyProduct(
                    jacobianInverseAdaptor.lazyProduct(
                        derivativeAdaptor));
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


template <maths::Expression TAD, maths::SpatialTransform TT>
unsigned LinearIsotropicStiffnessIntegrand<TAD,TT>::size() const noexcept
requires (!maths::StaticExpression<TAD>) {
    const auto derivativeComponentCount = _ansatzDerivatives.size();
    const auto ansatzCount = derivativeComponentCount / Dimension;
    return ansatzCount * ansatzCount;
}


template <maths::Expression TAD, maths::SpatialTransform TT>
constexpr unsigned LinearIsotropicStiffnessIntegrand<TAD,TT>::size() noexcept
requires (maths::StaticExpression<TAD>) {
    const auto derivativeComponentCount = TAD::size();
    const auto ansatzCount = derivativeComponentCount / Dimension;
    return ansatzCount * ansatzCount;
}


template <maths::Expression TAD, maths::SpatialTransform TT>
constexpr unsigned LinearIsotropicStiffnessIntegrand<TAD,TT>::bufferSize() noexcept
requires (maths::StaticExpression<TAD>) {
    return TAD::bufferSize();
}


template <maths::Expression TAD, maths::SpatialTransform TT>
unsigned LinearIsotropicStiffnessIntegrand<TAD,TT>::bufferSize() const noexcept
requires (!maths::StaticExpression<TAD>) {
    return _ansatzDerivatives.size() + _ansatzDerivatives.bufferSize();
}


} // namespace cie::fem
