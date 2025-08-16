#pragma once

// --- FEM Includes ---
#include "packages/maths/inc/Expression.hpp"

// --- STL Includes ---
#include <span> // span


namespace cie::fem::maths {


template <Expression TAnsatzDerivatives>
class LinearIsotropicStiffnessIntegrand
    : public ExpressionTraits<typename TAnsatzDerivatives::Value>
{
public:
    static constexpr unsigned Dimension = TAnsatzDerivatives::Dimension;

    using typename ExpressionTraits<typename TAnsatzDerivatives::Value>::Value;

    using typename ExpressionTraits<Value>::ConstSpan;

    using typename ExpressionTraits<Value>::Span;

public:
    LinearIsotropicStiffnessIntegrand();

    LinearIsotropicStiffnessIntegrand(const Value modulus,
                                      Ref<TAnsatzDerivatives> pAnsatzDerivatives);

    LinearIsotropicStiffnessIntegrand(const Value modulus,
                                      Ref<TAnsatzDerivatives> pAnsatzDerivatives,
                                      std::span<Value> buffer);

    void evaluate(ConstSpan in, Span out) const;

    unsigned size() const;

    unsigned getMinBufferSize() const noexcept;

    void setBuffer(std::span<Value> buffer);

private:
    Value _modulus;

    Ptr<TAnsatzDerivatives> _pAnsatzDerivatives;

    std::span<Value> _buffer;
}; // class LinearIsotropicStiffnessIntegrand


} // namespace cie::fem::maths

#include "packages/maths/impl/LinearIsotropicStiffnessIntegrand_impl.hpp"
