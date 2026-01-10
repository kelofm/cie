#pragma once

// --- FEM Includes ---
#include "packages/maths/inc/Expression.hpp"

// --- STL Includes ---
#include <span> // span


namespace cie::fem{


/// @details Computes
///          @f[
///              N_i \left( \xi \right)
///          @f]
///          and
///          @f[
///              N_i \left( \xi \right) u_d \left( x(\xi) \right)
///          @f]
///          where
///          - @f$ u_d @f$ is the desired state at global position @f$ x @f$,
///          - @f$ N_i @f$ is the ansatz function related to degree of freedom @f$ i @f$,
///          - @f$ \xi @f$ is the local coordinate within the cell to evaluate at,
///          - @f$ x @f$ is the global coordinate within the cell to evaluate at.
template <maths::Expression TDirichlet,
          maths::Expression TAnsatzSpace,
          maths::Expression TTransform>
class DirichletIntegrand
    : public maths::ExpressionTraits<typename TAnsatzSpace::Value> {
public:
    static constexpr unsigned Dimension = TAnsatzSpace::Dimension;

    using typename maths::ExpressionTraits<typename TAnsatzSpace::Value>::Value;

    using typename maths::ExpressionTraits<Value>::ConstSpan;

    using typename maths::ExpressionTraits<Value>::Span;

public:
    DirichletIntegrand();

    DirichletIntegrand(Ref<const TDirichlet> rDirichletFunctor,
                       Ref<const TAnsatzSpace> rAnsatzSpace,
                       Ref<const TTransform> rSpatialTransform);

    DirichletIntegrand(Ref<const TDirichlet> rDirichletFunctor,
                       Ref<const TAnsatzSpace> rAnsatzSpace,
                       Ref<const TTransform> rSpatialTransform,
                       std::span<Value> buffer);

    void evaluate(ConstSpan in, Span out) const;

    unsigned size() const;

    unsigned getMinBufferSize() const noexcept;

    void setBuffer(std::span<Value> buffer);

private:
    Ptr<const TDirichlet> _pDirichletFunctor;

    Ptr<const TAnsatzSpace> _pAnsatzSpace;

    Ptr<const TTransform> _pSpatialTransform;

    std::span<Value> _buffer;
}; // class DirichletIntegrand


template <maths::Expression TDirichlet,
          maths::Expression TAnsatzSpace,
          maths::Expression TTransform,
          concepts::Numeric TValue>
DirichletIntegrand<TDirichlet,TAnsatzSpace,TTransform>
makeDirichletIntegrand(Ref<const TDirichlet> rDirichletFunctor,
                       const TValue penalty,
                       Ref<const TAnsatzSpace> rAnsatzSpace,
                       Ref<const TTransform> rSpatialTransform,
                       std::span<TValue> buffer) {
    return DirichletIntegrand<TDirichlet,TAnsatzSpace,TTransform>(
        rDirichletFunctor,
        penalty,
        rAnsatzSpace,
        rSpatialTransform,
        buffer
    );
}


} // namespace cie::fem

#include "packages/integrands/impl/DirichletIntegrand_impl.hpp"
