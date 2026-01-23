#pragma once

// --- FEM Includes ---
#include "packages/utilities/inc/kernel.hpp"
#include "packages/maths/inc/Expression.hpp"

// --- Utility Includes ---
#include "packages/compile_time/packages/concepts/inc/basic_concepts.hpp"
#include "packages/stl_extension/inc/StaticArray.hpp"

// --- STL Includes ---
#include <algorithm>


namespace cie::fem {


template <unsigned Dim, concepts::Numeric TValue>
class QuadraturePoint {
public:
    using Value = typename Kernel<Dim,TValue>::LocalCoordinate;

    constexpr static inline unsigned Dimension = Dim;

    constexpr QuadraturePoint() noexcept {
        std::fill_n(
            _data.data(),
            Dim + 1,
            static_cast<TValue>(0));
    }

    constexpr QuadraturePoint(Ref<const std::span<const TValue,Dim>> rPosition,
                              TValue weight) noexcept {
        std::copy_n(
            rPosition.data(),
            Dim,
            _data.data());
        _data.back() = weight;
    }

    constexpr QuadraturePoint(Value position, Value weight) noexcept
    requires (Dim == 1u) {
        _data.front() = position;
        _data.back() = weight;
    }

    template <maths::Expression TExpression>
    requires std::is_same_v<typename TExpression::Value,TValue>
    void evaluate(Ref<const TExpression> rExpression,
                  std::span<TValue> out) const noexcept {
        rExpression.evaluate(Kernel<Dimension,TValue>::decay(this->position()), out);
        std::transform(
            out.begin(),
            out.end(),
            out.begin(),
            [weight = this->weight()] (TValue value) -> TValue {
                return value * weight;
            });
    }

    constexpr std::span<const Value,Dim> position() const noexcept {
        return std::span<const Value,Dim>(_data.data(), Dim);
    }

    constexpr Ref<const TValue> weight() const noexcept {
        return _data.back();
    }

    constexpr std::span<Value,Dim> position() noexcept {
        return std::span<Value,Dim>(_data.data(), Dim);
    }

    constexpr Ref<TValue> weight() noexcept {
        return _data.back();
    }

private:
    StaticArray<Value,Dim+1> _data;
};


} // namespace cie::fem
