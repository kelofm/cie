#pragma once

// --- FEM Includes ---
#include "packages/maths/inc/Expression.hpp"

// --- STL Includes ---
#include <span> // span


namespace cie::fem{


template <maths::Expression TAnsatzSpace>
class LinearIsotropicMassIntegrand
    : public maths::ExpressionTraits<typename TAnsatzSpace::Value>
{
public:
    static constexpr unsigned Dimension = TAnsatzSpace::Dimension;

    using typename maths::ExpressionTraits<typename TAnsatzSpace::Value>::Value;

    using typename maths::ExpressionTraits<Value>::ConstSpan;

    using typename maths::ExpressionTraits<Value>::Span;

public:
    LinearIsotropicMassIntegrand();

    LinearIsotropicMassIntegrand(const Value modulus,
                                 Ref<const TAnsatzSpace> rAnsatzSpace);

    LinearIsotropicMassIntegrand(const Value modulus,
                                 Ref<const TAnsatzSpace> rAnsatzSpace,
                                 std::span<Value> buffer);

    void evaluate(ConstSpan in, Span out) const;

    unsigned size() const;

    unsigned getMinBufferSize() const noexcept;

    void setBuffer(std::span<Value> buffer);

private:
    Value _modulus;

    Ptr<const TAnsatzSpace> _pAnsatzSpace;

    std::span<Value> _buffer;
}; // class LinearIsotropicMassIntegrand


} // namespace cie::fem

#include "packages/integrands/impl/LinearIsotropicMassIntegrand_impl.hpp"
