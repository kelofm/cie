#pragma once

// --- FEM Includes ---
#include "packages/maths/inc/Expression.hpp"


namespace cie::fem {


template <maths::Expression TAnsatzDerivatives, maths::SpatialTransform TTransform>
class LinearIsotropicStiffnessIntegrand
    : public maths::ExpressionTraits<typename TAnsatzDerivatives::Value> {
public:
    static constexpr unsigned ParametricDimension = TAnsatzDerivatives::ParametricDimension;

    using typename maths::ExpressionTraits<typename TAnsatzDerivatives::Value>::Value;

    using typename maths::ExpressionTraits<Value>::ConstSpan;

    using typename maths::ExpressionTraits<Value>::Span;

    using typename maths::ExpressionTraits<Value>::BufferSpan;

    using Jacobian = typename TTransform::Derivative;

    using JacobianInverse = typename TTransform::Inverse::Derivative;

    using Context = typename TAnsatzDerivatives::Context;

    using InContext = std::array<Context,2>;

    using OutContext = HyperSlab<2 * ParametricDimension,typename Context::Value>;

    template <template <class ...> class TOtherAllocator, class TOther>
    using Rebind = LinearIsotropicStiffnessIntegrand<
        typename TAnsatzDerivatives::template Rebind<TOtherAllocator,TOther>,
        typename TTransform::template Rebind<TOtherAllocator,TOther>>;

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

    /// @name Common Expression Interface
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

    void evaluate(
        ConstSpan in,
        Span out,
        BufferSpan buffer,
        Ref<const InContext> rInContext) const;

    /// @}
    /// @name Runtime Expression Interface
    /// @{

    unsigned size() const noexcept
    requires (!maths::StaticExpression<TAnsatzDerivatives>);

    unsigned bufferSize() const noexcept
    requires (!maths::StaticExpression<TAnsatzDerivatives>);

    /// @see @ref SliceableExpression
    unsigned size(Ref<const Context> rContext) const noexcept
    requires (!maths::StaticExpression<TAnsatzDerivatives>);

    unsigned size(Ref<const InContext> rInContext) const noexcept
    requires (!maths::StaticExpression<TAnsatzDerivatives>);

    /// @see @ref SliceableExpression
    unsigned bufferSize(Ref<const Context> rContext) const noexcept
    requires (!maths::StaticExpression<TAnsatzDerivatives>);

    unsigned bufferSize(Ref<const InContext> rInContext) const noexcept
    requires (!maths::StaticExpression<TAnsatzDerivatives>);

    /// @see @ref SliceableExpression
    std::array<typename Context::Value,ParametricDimension> contextBounds() const noexcept
    requires (!maths::StaticExpression<TAnsatzDerivatives>);

    OutContext outContext(Ref<const InContext> rInContext) const noexcept
    requires (!maths::StaticExpression<TAnsatzDerivatives>);

    /// @}
    /// @name Static Expression Interface
    /// @{

    static constexpr unsigned size() noexcept
    requires (maths::StaticExpression<TAnsatzDerivatives>);

    static constexpr unsigned bufferSize() noexcept
    requires (maths::StaticExpression<TAnsatzDerivatives>);

    /// @see @ref SliceableExpression
    static constexpr unsigned size(Ref<const Context> rContext) noexcept
    requires (maths::StaticExpression<TAnsatzDerivatives>);

    static constexpr unsigned size(Ref<const InContext> rInContext) noexcept
    requires (maths::StaticExpression<TAnsatzDerivatives>);

    /// @see @ref SliceableExpression
    static constexpr unsigned bufferSize(Ref<const Context> rContext) noexcept
    requires (maths::StaticExpression<TAnsatzDerivatives>);

    static constexpr unsigned bufferSize(Ref<const InContext> rInContext) noexcept
    requires (maths::StaticExpression<TAnsatzDerivatives>);

    /// @see @ref SliceableExpression
    static constexpr std::array<typename Context::Value,ParametricDimension> contextBounds() noexcept
    requires (maths::StaticExpression<TAnsatzDerivatives>);

    static constexpr OutContext outContext(Ref<const InContext> rInContext) noexcept
    requires (maths::StaticExpression<TAnsatzDerivatives>);

    /// @}
    /// @name Serialization
    /// @{

    void serialize(
        Ref<cie::io::Traits::SerializerStream> rStream,
        tags::Binary tag = {}) const;

    template <class TAllocator>
    static void deserialize(
        Ref<cie::io::Traits::DeserializerStream> rStream,
        Ref<LinearIsotropicStiffnessIntegrand> rInstance,
        Ref<const TAllocator> rAllocator,
        tags::Binary tag = {});

    /// @}

private:
    TAnsatzDerivatives _ansatzDerivatives;

    Jacobian _jacobian;

    JacobianInverse _jacobianInverse;

    Value _modulus;
}; // class LinearIsotropicStiffnessIntegrand


} // namespace cie::fem

#include "packages/integrands/impl/LinearIsotropicStiffnessIntegrand_impl.hpp"
