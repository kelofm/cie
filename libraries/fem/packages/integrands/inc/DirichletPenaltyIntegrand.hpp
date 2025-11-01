#pragma once

// --- FEM Includes ---
#include "packages/maths/inc/Expression.hpp"

// --- STL Includes ---
#include <span> // span


namespace cie::fem{


template <maths::Expression TDirichlet, maths::Expression TAnsatzSpace>
class DirichletPenaltyIntegrand
    : public maths::ExpressionTraits<typename TAnsatzSpace::Value>
{
public:
    static constexpr unsigned Dimension = TAnsatzSpace::Dimension;

    using typename maths::ExpressionTraits<typename TAnsatzSpace::Value>::Value;

    using typename maths::ExpressionTraits<Value>::ConstSpan;

    using typename maths::ExpressionTraits<Value>::Span;

public:
    DirichletPenaltyIntegrand();

    DirichletPenaltyIntegrand(Ref<const TDirichlet> rDirichletFunctor,
                              const Value penalty,
                              Ref<const TAnsatzSpace> rAnsatzSpace);

    DirichletPenaltyIntegrand(Ref<const TDirichlet> rDirichletFunctor,
                              const Value penalty,
                              Ref<const TAnsatzSpace> rAnsatzSpace,
                              std::span<Value> buffer);

    void evaluate(ConstSpan in, Span out) const;

    unsigned size() const;

    unsigned getMinBufferSize() const noexcept;

    void setBuffer(std::span<Value> buffer);

private:
    Value _penalty;

    Ptr<const TDirichlet> _pDirichletFunctor;

    Ptr<const TAnsatzSpace> _pAnsatzSpace;

    std::span<Value> _buffer;
}; // class DirichletPenaltyIntegrand


template <maths::Expression TDirichlet,
          maths::Expression TAnsatzSpace,
          concepts::Numeric TValue>
DirichletPenaltyIntegrand<TDirichlet,TAnsatzSpace>
makeDirichletPenaltyIntegrand(Ref<const TDirichlet> rDirichletFunctor,
                              const TValue penalty,
                              Ref<const TAnsatzSpace> rAnsatzSpace,
                              std::span<TValue> buffer)
{
    return DirichletPenaltyIntegrand<TDirichlet,TAnsatzSpace>(
        rDirichletFunctor,
        penalty,
        rAnsatzSpace,
        buffer
    );
}


} // namespace cie::fem

#include "packages/integrands/impl/DirichletPenaltyIntegrand_impl.hpp"
