#pragma once

// --- FEM Includes ---
#include "packages/maths/inc/Expression.hpp"
#include "packages/numeric/inc/CellBase.hpp"

// --- Utility Includes ---
#include "packages/io/inc/Traits.hpp"

// --- STL Includes ---
#include <array>


namespace cie::fem{


/// @brief Legacy combined Dirichlet penalty contribution.
/// @details Kept as a reference implementation; prefer @ref DirichletPenaltyLHSIntegrand and
///          @ref DirichletPenaltyRHSIntegrand when assembling sliced contributions.
///          Computes
///          @f[
///              N_i \left( \xi \right) p N_j \left(\xi \right)
///          @f]
///          and
///          @f[
///              p N_i \left( \xi \right) u_d \left( x(\xi) \right)
///          @f]
///          if @p Symmetric is true, otherwise
///          @f[
///              p N_i \left( \xi \right)
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
template <
    maths::Expression TDirichlet,
    maths::SliceableExpression TAnsatzSpace,
    maths::Expression TEmbedding,
    CellLike TCell>
class DirichletPenaltyIntegrand : public maths::ExpressionTraits<typename TAnsatzSpace::Value> {
public:
    static constexpr unsigned ParametricDimension = TAnsatzSpace::ParametricDimension;

    using typename maths::ExpressionTraits<typename TAnsatzSpace::Value>::Value;

    using typename maths::ExpressionTraits<Value>::ConstSpan;

    using typename maths::ExpressionTraits<Value>::Span;

    using typename maths::ExpressionTraits<Value>::BufferSpan;

    using CellInverseTransform = typename TCell::SpatialTransform::Inverse;

    using CellJacobian = typename TCell::SpatialTransform::Derivative;

    using Context = TAnsatzSpace::Context;

    static inline constexpr bool isStatic = (
           maths::StaticExpression<TDirichlet>
        && maths::StaticExpression<TAnsatzSpace>
        && maths::StaticExpression<TEmbedding>
        && maths::StaticExpression<CellInverseTransform>);

    template <template <class ...> class TOtherAllocator, class TOther>
    using Rebind = DirichletPenaltyIntegrand<
        typename TDirichlet::template Rebind<TOtherAllocator,TOther>,
        typename TAnsatzSpace::template Rebind<TOtherAllocator,TOther>,
        typename TEmbedding::template Rebind<TOtherAllocator,TOther>,
        TCell>;

public:
    DirichletPenaltyIntegrand() noexcept;

    DirichletPenaltyIntegrand(
        Ref<const TDirichlet> rDirichletFunctor,
        const Value penalty,
        Ref<const TAnsatzSpace> rAnsatzSpace,
        Ref<const TEmbedding> rSpatialTransform,
        Ref<const CellInverseTransform> rCellInverseTransform,
        Ref<const CellJacobian> rCellJacobian);

    /// @name Runtime Expression Interface
    /// @{

    void evaluate(
        ConstSpan in,
        Span out,
        BufferSpan buffer) const;

    void evaluate(
        ConstSpan in,
        Span out,
        BufferSpan buffer,
        Ref<const Context> rContext) const;

    unsigned size() const noexcept;

    unsigned bufferSize() const noexcept;

    unsigned size(Ref<const Context> rContext) const noexcept;

    unsigned bufferSize(Ref<const Context> rContext) const noexcept;

    std::array<typename Context::Value,ParametricDimension> contextBounds() const noexcept;

    /// @}
    /// @name Serialization
    /// @{

    void serialize(
        Ref<cie::io::Traits::SerializerStream> rStream,
        tags::Binary tag = {}) const;

    template <class TStream, class TAllocator>
    static void deserialize(
        Ref<TStream> rStream,
        Ref<DirichletPenaltyIntegrand> rInstance,
        Ref<const TAllocator> rAllocator,
        tags::Binary tag = {});

    /// @}

private:
    TAnsatzSpace _ansatzSpace;

    TEmbedding _embedding;

    CellInverseTransform _cellInverseTransform;

    CellJacobian _cellJacobian;

    TDirichlet _dirichletFunctor;

    Value _penalty;
}; // class DirichletPenaltyIntegrand


/// @brief Dirichlet penalty contribution to the local left-hand-side matrix.
/// @details The first input context selects rows and the second selects columns.
template <
    maths::SliceableExpression TAnsatzSpace,
    maths::Expression TEmbedding,
    CellLike TCell>
class DirichletPenaltyLHSIntegrand : public maths::ExpressionTraits<typename TAnsatzSpace::Value> {
public:
    static constexpr unsigned ParametricDimension = TAnsatzSpace::ParametricDimension;

    using typename maths::ExpressionTraits<typename TAnsatzSpace::Value>::Value;
    using typename maths::ExpressionTraits<Value>::ConstSpan;
    using typename maths::ExpressionTraits<Value>::Span;
    using typename maths::ExpressionTraits<Value>::BufferSpan;

    using Context = typename TAnsatzSpace::Context;
    using InContext = std::array<Context,2>;
    using OutContext = HyperSlab<2 * ParametricDimension,typename Context::Value>;
    using CellInverseTransform = typename TCell::SpatialTransform::Inverse;
    using CellJacobian = typename TCell::SpatialTransform::Derivative;

    template <template <class ...> class TOtherAllocator, class TOther>
    using Rebind = DirichletPenaltyLHSIntegrand<
        typename TAnsatzSpace::template Rebind<TOtherAllocator,TOther>,
        typename TEmbedding::template Rebind<TOtherAllocator,TOther>,
        TCell>;

    DirichletPenaltyLHSIntegrand() noexcept;

    DirichletPenaltyLHSIntegrand(
        const Value penalty,
        Ref<const TAnsatzSpace> rAnsatzSpace,
        Ref<const TEmbedding> rSpatialTransform,
        Ref<const CellInverseTransform> rCellInverseTransform,
        Ref<const CellJacobian> rCellJacobian);

    void evaluate(ConstSpan in, Span out, BufferSpan buffer) const;

    void evaluate(
        ConstSpan in,
        Span out,
        BufferSpan buffer,
        Ref<const InContext> rInContext) const;

    unsigned size() const noexcept;
    unsigned size(Ref<const InContext> rInContext) const noexcept;
    unsigned bufferSize() const noexcept;
    unsigned bufferSize(Ref<const InContext> rInContext) const noexcept;
    std::array<typename Context::Value,ParametricDimension> contextBounds() const noexcept;
    OutContext outContext(Ref<const InContext> rInContext) const noexcept;

    void serialize(
        Ref<cie::io::Traits::SerializerStream> rStream,
        tags::Binary tag = {}) const;

    template <class TStream, class TAllocator>
    static void deserialize(
        Ref<TStream> rStream,
        Ref<DirichletPenaltyLHSIntegrand> rInstance,
        Ref<const TAllocator> rAllocator,
        tags::Binary tag = {});

private:
    TAnsatzSpace _ansatzSpace;
    TEmbedding _embedding;
    CellInverseTransform _cellInverseTransform;
    CellJacobian _cellJacobian;
    Value _penalty;
}; // class DirichletPenaltyLHSIntegrand


/// @brief Dirichlet penalty contribution to the local right-hand-side vector.
template <
    maths::Expression TDirichlet,
    maths::SliceableExpression TAnsatzSpace,
    maths::Expression TEmbedding,
    CellLike TCell>
class DirichletPenaltyRHSIntegrand : public maths::ExpressionTraits<typename TAnsatzSpace::Value> {
public:
    static constexpr unsigned ParametricDimension = TAnsatzSpace::ParametricDimension;

    using typename maths::ExpressionTraits<typename TAnsatzSpace::Value>::Value;
    using typename maths::ExpressionTraits<Value>::ConstSpan;
    using typename maths::ExpressionTraits<Value>::Span;
    using typename maths::ExpressionTraits<Value>::BufferSpan;

    using Context = typename TAnsatzSpace::Context;
    using InContext = std::array<Context,1>;
    using OutContext = HyperSlab<ParametricDimension+1,typename Context::Value>;
    using CellInverseTransform = typename TCell::SpatialTransform::Inverse;
    using CellJacobian = typename TCell::SpatialTransform::Derivative;

    template <template <class ...> class TOtherAllocator, class TOther>
    using Rebind = DirichletPenaltyRHSIntegrand<
        typename TDirichlet::template Rebind<TOtherAllocator,TOther>,
        typename TAnsatzSpace::template Rebind<TOtherAllocator,TOther>,
        typename TEmbedding::template Rebind<TOtherAllocator,TOther>,
        TCell>;

    DirichletPenaltyRHSIntegrand() noexcept;

    DirichletPenaltyRHSIntegrand(
        Ref<const TDirichlet> rDirichletFunctor,
        const Value penalty,
        Ref<const TAnsatzSpace> rAnsatzSpace,
        Ref<const TEmbedding> rSpatialTransform,
        Ref<const CellInverseTransform> rCellInverseTransform,
        Ref<const CellJacobian> rCellJacobian);

    void evaluate(ConstSpan in, Span out, BufferSpan buffer) const;

    void evaluate(
        ConstSpan in,
        Span out,
        BufferSpan buffer,
        Ref<const InContext> rInContext) const;

    unsigned size() const noexcept;
    unsigned size(Ref<const InContext> rInContext) const noexcept;
    unsigned bufferSize() const noexcept;
    unsigned bufferSize(Ref<const InContext> rInContext) const noexcept;
    std::array<typename Context::Value,ParametricDimension> contextBounds() const noexcept;
    OutContext outContext(Ref<const InContext> rInContext) const noexcept;

    void serialize(
        Ref<cie::io::Traits::SerializerStream> rStream,
        tags::Binary tag = {}) const;

    template <class TStream, class TAllocator>
    static void deserialize(
        Ref<TStream> rStream,
        Ref<DirichletPenaltyRHSIntegrand> rInstance,
        Ref<const TAllocator> rAllocator,
        tags::Binary tag = {});

private:
    TAnsatzSpace _ansatzSpace;
    TEmbedding _embedding;
    CellInverseTransform _cellInverseTransform;
    CellJacobian _cellJacobian;
    TDirichlet _dirichletFunctor;
    Value _penalty;
}; // class DirichletPenaltyRHSIntegrand


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
    Ref<const TCell> rCell) {
        return DirichletPenaltyIntegrand<TDirichlet,TAnsatzSpace,TEmbedding,TCell>(
            rDirichletFunctor,
            penalty,
            rAnsatzSpace,
            rSpatialTransform,
            rCell.makeInverseSpatialTransform(),
            rCell.makeJacobian());
}


template <maths::SliceableExpression TAnsatzSpace,
          maths::Expression TEmbedding,
          CellLike TCell,
          concepts::Numeric TValue>
DirichletPenaltyLHSIntegrand<TAnsatzSpace,TEmbedding,TCell>
makeDirichletPenaltyLHSIntegrand(
    const TValue penalty,
    Ref<const TAnsatzSpace> rAnsatzSpace,
    Ref<const TEmbedding> rSpatialTransform,
    Ref<const TCell> rCell) {
        return DirichletPenaltyLHSIntegrand<TAnsatzSpace,TEmbedding,TCell>(
            penalty,
            rAnsatzSpace,
            rSpatialTransform,
            rCell.makeInverseSpatialTransform(),
            rCell.makeJacobian());
}


template <maths::Expression TDirichlet,
          maths::SliceableExpression TAnsatzSpace,
          maths::Expression TEmbedding,
          CellLike TCell,
          concepts::Numeric TValue>
DirichletPenaltyRHSIntegrand<TDirichlet,TAnsatzSpace,TEmbedding,TCell>
makeDirichletPenaltyRHSIntegrand(
    Ref<const TDirichlet> rDirichletFunctor,
    const TValue penalty,
    Ref<const TAnsatzSpace> rAnsatzSpace,
    Ref<const TEmbedding> rSpatialTransform,
    Ref<const TCell> rCell) {
        return DirichletPenaltyRHSIntegrand<TDirichlet,TAnsatzSpace,TEmbedding,TCell>(
            rDirichletFunctor,
            penalty,
            rAnsatzSpace,
            rSpatialTransform,
            rCell.makeInverseSpatialTransform(),
            rCell.makeJacobian());
}


} // namespace cie::fem

#include "packages/integrands/impl/DirichletPenaltyIntegrand_impl.hpp"
