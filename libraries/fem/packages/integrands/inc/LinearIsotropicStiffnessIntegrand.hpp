#pragma once

// --- FEM Includes ---
#include "packages/maths/inc/Expression.hpp"

// --- STL Includes ---
#include <span> // span


namespace cie::fem {


template <maths::Expression TAnsatzDerivatives>
class LinearIsotropicStiffnessIntegrand
    : public maths::ExpressionTraits<typename TAnsatzDerivatives::Value> {
public:
    static constexpr unsigned Dimension = TAnsatzDerivatives::Dimension;

    using typename maths::ExpressionTraits<typename TAnsatzDerivatives::Value>::Value;

    using typename maths::ExpressionTraits<Value>::ConstSpan;

    using typename maths::ExpressionTraits<Value>::Span;

public:
    static inline constexpr bool isBuffered = maths::BufferedExpression<TAnsatzDerivatives> || !maths::StaticExpression<TAnsatzDerivatives>;

    LinearIsotropicStiffnessIntegrand();

    LinearIsotropicStiffnessIntegrand(
        const Value modulus,
        RightRef<TAnsatzDerivatives> rAnsatzDerivatives) noexcept;

    LinearIsotropicStiffnessIntegrand(
        const Value modulus,
        RightRef<TAnsatzDerivatives> rAnsatzDerivatives,
        std::span<Value> buffer)
    requires (isBuffered);

    LinearIsotropicStiffnessIntegrand(
        const Value modulus,
        Ref<const TAnsatzDerivatives> rAnsatzDerivatives);

    LinearIsotropicStiffnessIntegrand(
        const Value modulus,
        Ref<const TAnsatzDerivatives> rAnsatzDerivatives,
        std::span<Value> buffer)
    requires (isBuffered);

    void evaluate(ConstSpan in, Span out) const;

    unsigned size() const noexcept
    requires (!maths::StaticExpression<TAnsatzDerivatives>);

    static constexpr unsigned size() noexcept
    requires (maths::StaticExpression<TAnsatzDerivatives>);

    unsigned getMinBufferSize() const noexcept
    requires (isBuffered);

    void setBuffer(std::span<Value> buffer)
    requires (isBuffered);

private:
    Value _modulus;

    TAnsatzDerivatives _ansatzDerivatives;

    mutable std::conditional_t<
        isBuffered,
        std::span<Value>,
        std::array<Value,maths::StaticExpressionSize<TAnsatzDerivatives>::value>
    > _buffer;
}; // class LinearIsotropicStiffnessIntegrand


} // namespace cie::fem

#include "packages/integrands/impl/LinearIsotropicStiffnessIntegrand_impl.hpp"
