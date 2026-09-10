#pragma once

// --- FEM Includes ---
#include "packages/maths/inc/Expression.hpp"

// --- STL Includes ---
#include <array>
#include <utility>


namespace cie::fem {


/// @ingroup fem
template <maths::SliceableExpression TIntegrand, class TMaterialID, int MaterialCount>
class ScaledMultiMaterialIntegrand : public maths::ExpressionTraits<typename TIntegrand::Value> {
public:
    static constexpr unsigned ParametricDimension = TIntegrand::ParametricDimension;

    using typename maths::ExpressionTraits<typename TIntegrand::Value>::Value;

    using typename maths::ExpressionTraits<Value>::ConstSpan;

    using typename maths::ExpressionTraits<Value>::Span;

    using typename maths::ExpressionTraits<Value>::BufferSpan;

    using Context = typename TIntegrand::Context;

    template <template <class ...> class TOtherAllocator, class TOther>
    using Rebind = ScaledMultiMaterialIntegrand<
        typename TIntegrand::template Rebind<TOtherAllocator,TOther>,
        TMaterialID,
        MaterialCount>;

    ScaledMultiMaterialIntegrand() noexcept = default;

    ScaledMultiMaterialIntegrand(
        RightRef<TIntegrand> rIntegrand,
        std::span<const std::pair<TMaterialID,Value>,MaterialCount> materialMap);

    /// @name Common Expression Interface
    /// @{

    void evaluate(
        ConstSpan in,
        Span out,
        BufferSpan buffer) const;

    /// @see @ref SliceableExpression
    void evaluate(
        ConstSpan in,
        Span out,
        BufferSpan buffer,
        Ref<const Context> rContext) const;

    /// @}
    /// @name Runtime Expression Interface
    /// @{

    unsigned size() const noexcept
    requires (!maths::StaticExpression<TIntegrand>);

    unsigned bufferSize() const noexcept
    requires (!maths::StaticExpression<TIntegrand>);

    /// @see @ref SliceableExpression
    unsigned size(Ref<const Context> rContext) const noexcept
    requires (!maths::StaticExpression<TIntegrand>);

    /// @see @ref SliceableExpression
    unsigned bufferSize(Ref<const Context> rContext)  const noexcept
    requires (!maths::StaticExpression<TIntegrand>);

    /// @see @ref SliceableExpression
    std::array<typename Context::Value,ParametricDimension> contextBounds() const noexcept
    requires (!maths::StaticExpression<TIntegrand>);

    /// @}
    /// @name Static Expression Interface
    /// @{

    static constexpr unsigned size() noexcept
    requires maths::StaticExpression<TIntegrand>;

    static constexpr unsigned bufferSize() noexcept
    requires maths::StaticExpression<TIntegrand>;

    /// @see @ref SliceableExpression
    static constexpr unsigned size(Ref<const Context> rContext) noexcept
    requires maths::StaticExpression<TIntegrand>;

    /// @see @ref SliceableExpression
    static constexpr unsigned bufferSize(Ref<const Context> rContext) noexcept
    requires maths::StaticExpression<TIntegrand>;

    /// @see @ref SliceableExpression
    static constexpr std::array<typename Context::Value,ParametricDimension> contextBounds() noexcept
    requires (maths::StaticExpression<TIntegrand>);

    /// @}
    /// @name Serialization
    /// @{

    void serialize(
        Ref<cie::io::Traits::SerializerStream> rStream,
        tags::Binary tag = {}) const;

    static void deserialize(
        Ref<cie::io::Traits::DeserializerStream> rStream,
        Ref<ScaledMultiMaterialIntegrand> rInstance,
        tags::Binary tag = {});

    /// @}

private:
    TIntegrand _integrand;

    std::array<
        std::pair<
            TMaterialID,
            Value
        >,
        MaterialCount
    > _materialMap;
}; // class ScaledMultiMaterialIntegrand


} // namespace cie::fem

#include "packages/integrands/impl/ScaledMultiMaterialIntegrand_impl.hpp"
