#pragma once

// --- FEM Includes ---
#include "packages/maths/inc/Expression.hpp"
#include "packages/numeric/inc/Cell.hpp"

// --- STL Includes ---
#include <span> // span


namespace cie::fem{


/// @details Computes
///          @f[
///              N_i \left( \xi \right) p N_j \left(\xi \right)
///          @f]
///          and
///          @f[
///              p N_i \left( \xi \right) u_d \left( x(\xi) \right)
///          @f]
///          where
///          - @f$ p @f$ is the penalty factor,
///          - @f$ u_d @f$ is the desired state at global position @f$ x @f$,
///          - @f$ N_i @f$ is the ansatz function related to degree of freedom @f$ i @f$,
///          - @f$ \xi @f$ is the local coordinate within the cell to evaluate at,
///          - @f$ x @f$ is the global coordinate within the cell to evaluate at.
template <maths::Expression TDirichlet,
          maths::Expression TAnsatzSpace,
          maths::Expression TEmbedding,
          CellLike TCell>
class DirichletPenaltyIntegrand
    : public maths::ExpressionTraits<typename TAnsatzSpace::Value>
{
public:
    static constexpr unsigned Dimension = TAnsatzSpace::Dimension;

    using typename maths::ExpressionTraits<typename TAnsatzSpace::Value>::Value;

    using typename maths::ExpressionTraits<Value>::ConstSpan;

    using typename maths::ExpressionTraits<Value>::Span;

    using CellInverseTransform = typename TCell::SpatialTransform::Inverse;

    static inline constexpr bool isBuffered = (
           maths::BufferedExpression<TDirichlet>
        || maths::BufferedExpression<TAnsatzSpace>
        || maths::BufferedExpression<TEmbedding>
        || maths::BufferedExpression<CellInverseTransform>);

    static inline constexpr bool isStatic = (
           maths::StaticExpression<TDirichlet>
        || maths::StaticExpression<TAnsatzSpace>
        || maths::StaticExpression<TEmbedding>
        || maths::StaticExpression<CellInverseTransform>);

public:
    DirichletPenaltyIntegrand();

    DirichletPenaltyIntegrand(
        Ref<const TDirichlet> rDirichletFunctor,
        const Value penalty,
        Ref<const TAnsatzSpace> rAnsatzSpace,
        Ref<const TEmbedding> rSpatialTransform,
        Ref<const CellInverseTransform> rCellInverseTransform);

    DirichletPenaltyIntegrand(
        Ref<const TDirichlet> rDirichletFunctor,
        const Value penalty,
        Ref<const TAnsatzSpace> rAnsatzSpace,
        Ref<const TEmbedding> rSpatialTransform,
        Ref<const CellInverseTransform> rCellInverseTransform,
        std::span<Value> buffer);

    void evaluate(ConstSpan in, Span out) const;

    unsigned size() const;

    unsigned getMinBufferSize() const noexcept;

    void setBuffer(std::span<Value> buffer);

private:
    Value _penalty;

    Ptr<const TDirichlet> _pDirichletFunctor;

    Ptr<const TAnsatzSpace> _pAnsatzSpace;

    Ptr<const TEmbedding> _pEmbedding;

    Ptr<const CellInverseTransform> _pCellInverseTransform;

    std::span<Value> _buffer;
}; // class DirichletPenaltyIntegrand


template <maths::Expression TDirichlet,
          maths::Expression TAnsatzSpace,
          maths::Expression TEmbedding,
          CellLike TCell,
          concepts::Numeric TValue>
DirichletPenaltyIntegrand<TDirichlet,TAnsatzSpace,TEmbedding,TCell>
makeDirichletPenaltyIntegrand(
    Ref<const TDirichlet> rDirichletFunctor,
    const TValue penalty,
    Ref<const TAnsatzSpace> rAnsatzSpace,
    Ref<const TEmbedding> rSpatialTransform,
    Ref<const TCell> rCell,
    std::span<TValue> buffer = {})
{
    if (buffer.empty()) {
        return DirichletPenaltyIntegrand<TDirichlet,TAnsatzSpace,TEmbedding,TCell>(
            rDirichletFunctor,
            penalty,
            rAnsatzSpace,
            rSpatialTransform,
            rCell.makeInverseSpatialTransform());
    } else {
        return DirichletPenaltyIntegrand<TDirichlet,TAnsatzSpace,TEmbedding,TCell>(
            rDirichletFunctor,
            penalty,
            rAnsatzSpace,
            rSpatialTransform,
            rCell.makeInverseSpatialTransform(),
            buffer);
    }
}


} // namespace cie::fem

#include "packages/integrands/impl/DirichletPenaltyIntegrand_impl.hpp"
