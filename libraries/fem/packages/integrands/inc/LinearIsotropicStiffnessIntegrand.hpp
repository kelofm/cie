#pragma once

// --- FEM Includes ---
#include "packages/maths/inc/Expression.hpp"


namespace cie::fem {


template <maths::Expression TAnsatzDerivatives>
class LinearIsotropicStiffnessIntegrand
    : public maths::ExpressionTraits<typename TAnsatzDerivatives::Value> {
public:
    static constexpr unsigned Dimension = TAnsatzDerivatives::Dimension;

    using typename maths::ExpressionTraits<typename TAnsatzDerivatives::Value>::Value;

    using typename maths::ExpressionTraits<Value>::ConstSpan;

    using typename maths::ExpressionTraits<Value>::Span;

    using typename maths::ExpressionTraits<Value>::BufferSpan;

public:
    LinearIsotropicStiffnessIntegrand();

    LinearIsotropicStiffnessIntegrand(
        const Value modulus,
        RightRef<TAnsatzDerivatives> rAnsatzDerivatives) noexcept;

    LinearIsotropicStiffnessIntegrand(
        const Value modulus,
        Ref<const TAnsatzDerivatives> rAnsatzDerivatives);

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

private:
    Value _modulus;

    TAnsatzDerivatives _ansatzDerivatives;
}; // class LinearIsotropicStiffnessIntegrand


} // namespace cie::fem

#include "packages/integrands/impl/LinearIsotropicStiffnessIntegrand_impl.hpp"
