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
    LinearIsotropicStiffnessIntegrand();

    LinearIsotropicStiffnessIntegrand(const Value modulus,
                                      RightRef<TAnsatzDerivatives> rAnsatzDerivatives)
    requires maths::BufferedExpression<TAnsatzDerivatives>;

    LinearIsotropicStiffnessIntegrand(const Value modulus,
                                      RightRef<TAnsatzDerivatives> rAnsatzDerivatives,
                                      std::span<Value> buffer)
    requires maths::BufferedExpression<TAnsatzDerivatives>;

    LinearIsotropicStiffnessIntegrand(const Value modulus,
                                      Ref<const TAnsatzDerivatives> rAnsatzDerivatives)
    requires (!maths::BufferedExpression<TAnsatzDerivatives>);

    LinearIsotropicStiffnessIntegrand(const Value modulus,
                                      Ref<const TAnsatzDerivatives> rAnsatzDerivatives,
                                      std::span<Value> buffer)
    requires (!maths::BufferedExpression<TAnsatzDerivatives>);

    void evaluate(ConstSpan in, Span out) const;

    unsigned size() const;

    unsigned getMinBufferSize() const noexcept;

    void setBuffer(std::span<Value> buffer);

private:
    Value _modulus;

    std::conditional_t<
        maths::BufferedExpression<TAnsatzDerivatives>,
        TAnsatzDerivatives,
        Ptr<const TAnsatzDerivatives>
    > _ansatzDerivatives;

    std::span<Value> _buffer;
}; // class LinearIsotropicStiffnessIntegrand


} // namespace cie::fem

#include "packages/integrands/impl/LinearIsotropicStiffnessIntegrand_impl.hpp"
