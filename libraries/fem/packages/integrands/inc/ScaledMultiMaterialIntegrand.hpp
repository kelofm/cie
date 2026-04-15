#pragma once

// --- FEM Includes ---
#include "packages/maths/inc/Expression.hpp"

// --- STL Includes ---
#include <array>
#include <utility>


namespace cie::fem {


template <maths::Expression TIntegrand, class TMaterialID, int MaterialCount>
class ScaledMultiMaterialIntegrand : public maths::ExpressionTraits<typename TIntegrand::Value> {
public:
    static constexpr unsigned Dimension = TIntegrand::Dimension;

    using typename maths::ExpressionTraits<typename TIntegrand::Value>::Value;

    using typename maths::ExpressionTraits<Value>::ConstSpan;

    using typename maths::ExpressionTraits<Value>::Span;

    using typename maths::ExpressionTraits<Value>::BufferSpan;

    ScaledMultiMaterialIntegrand() noexcept = default;

    ScaledMultiMaterialIntegrand(
        RightRef<TIntegrand> rIntegrand,
        std::span<const std::pair<TMaterialID,Value>,MaterialCount> materialMap);

    void evaluate(
        ConstSpan in,
        Span out,
        BufferSpan buffer) const;

    unsigned size() const noexcept
    requires (!maths::StaticExpression<TIntegrand>);

    static constexpr unsigned size() noexcept
    requires maths::StaticExpression<TIntegrand>;

    unsigned bufferSize() const noexcept
    requires (!maths::StaticExpression<TIntegrand>);

    static constexpr unsigned bufferSize() noexcept
    requires maths::StaticExpression<TIntegrand>;

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
