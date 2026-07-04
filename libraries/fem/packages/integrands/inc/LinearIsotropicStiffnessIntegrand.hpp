#pragma once

// --- FEM Includes ---
#include "packages/maths/inc/Expression.hpp"


namespace cie::fem {


template <maths::Expression TAnsatzDerivatives, maths::SpatialTransform TTransform>
class LinearIsotropicStiffnessIntegrand
    : public maths::ExpressionTraits<typename TAnsatzDerivatives::Value> {
public:
    static constexpr unsigned Dimension = TAnsatzDerivatives::Dimension;

    using typename maths::ExpressionTraits<typename TAnsatzDerivatives::Value>::Value;

    using typename maths::ExpressionTraits<Value>::ConstSpan;

    using typename maths::ExpressionTraits<Value>::Span;

    using typename maths::ExpressionTraits<Value>::BufferSpan;

    using Jacobian = typename TTransform::Derivative;

    using JacobianInverse = typename TTransform::Inverse::Derivative;

public:
    LinearIsotropicStiffnessIntegrand();

    LinearIsotropicStiffnessIntegrand(
        const Value modulus,
        RightRef<TAnsatzDerivatives> rAnsatzDerivatives,
        RightRef<Jacobian> rJacobian,
        RightRef<JacobianInverse> rJacobianInverse) noexcept;

    LinearIsotropicStiffnessIntegrand(
        const Value modulus,
        Ref<const TAnsatzDerivatives> rAnsatzDerivatives,
        Ref<const Jacobian> rJacobian,
        Ref<const JacobianInverse> rJacobianInverse);

    void evaluate(
        ConstSpan in,
        Span out,
        BufferSpan buffer) const;

    unsigned size() const noexcept
    requires (!maths::StaticExpression<TAnsatzDerivatives>);

    static constexpr unsigned size() noexcept
    requires (maths::StaticExpression<TAnsatzDerivatives>);

    unsigned bufferSize() const noexcept
    requires (!maths::StaticExpression<TAnsatzDerivatives>);

    static constexpr unsigned bufferSize() noexcept
    requires (maths::StaticExpression<TAnsatzDerivatives>);

    void serialize(
        Ref<cie::io::Traits::SerializerStream> rStream,
        tags::Binary tag = {}) const;

    static void deserialize(
        Ref<cie::io::Traits::DeserializerStream> rStream,
        Ref<LinearIsotropicStiffnessIntegrand> rInstance,
        tags::Binary tag = {});

private:
    TAnsatzDerivatives _ansatzDerivatives;

    Jacobian _jacobian;

    JacobianInverse _jacobianInverse;

    Value _modulus;
}; // class LinearIsotropicStiffnessIntegrand


} // namespace cie::fem

#include "packages/integrands/impl/LinearIsotropicStiffnessIntegrand_impl.hpp"
