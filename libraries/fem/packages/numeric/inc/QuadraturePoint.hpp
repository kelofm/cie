#pragma once

// --- FEM Includes ---
#include "packages/utilities/inc/kernel.hpp"
#include "packages/maths/inc/Expression.hpp"
#include "packages/maths/inc/IdentityTransform.hpp"

// --- Utility Includes ---
#include "packages/compile_time/packages/concepts/inc/basic_concepts.hpp"
#include "packages/stl_extension/inc/StaticArray.hpp"

// --- STL Includes ---
#include <algorithm>


namespace cie::fem {


template <
    unsigned Dim,
    concepts::Numeric TValue,
    class TData = void>
class QuadraturePoint {
public:
    using Value = ParametricCoordinate<TValue>;

    constexpr static inline unsigned Dimension = Dim;

    constexpr QuadraturePoint() noexcept {
        std::fill_n(
            std::get<0>(_data).data(),
            Dim + 1,
            static_cast<TValue>(0));
    }

    constexpr QuadraturePoint(Ref<const std::span<const TValue,Dim>> rPosition,
                              TValue weight) noexcept {
        std::copy_n(
            rPosition.data(),
            Dim,
            std::get<0>(_data).data());
        std::get<0>(_data).back() = weight;
    }

    constexpr QuadraturePoint(
        Ref<const std::span<const TValue,Dim>> rPosition,
        TValue weight,
        typename VoidSafe<TData,int>::RightRef rData) noexcept
    requires (!std::is_same_v<TData,void>)
        : QuadraturePoint(rPosition, weight) {
        this->data() = std::move(rData);
    }

    constexpr QuadraturePoint(Value position, Value weight) noexcept
    requires (Dim == 1u) {
        std::get<0>(_data).front() = position;
        std::get<0>(_data).back() = weight;
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
        return std::span<const Value,Dim>(std::get<0>(_data).data(), Dim);
    }

    constexpr Ref<const TValue> weight() const noexcept {
        return std::get<0>(_data).back();
    }

    constexpr std::span<Value,Dim> position() noexcept {
        return std::span<Value,Dim>(std::get<0>(_data).data(), Dim);
    }

    constexpr Ref<TValue> weight() noexcept {
        return std::get<0>(_data).back();
    }

    constexpr typename VoidSafe<const TData>::Ref data() const noexcept {
        if constexpr (std::is_same_v<TData,void>) return;
        else return std::get<1>(_data);
    }

    constexpr typename VoidSafe<TData>::Ref data() noexcept {
        if constexpr (std::is_same_v<TData,void>) return;
        else return std::get<1>(_data);
    }

private:
    std::conditional_t<
        std::is_same_v<TData,void>,
        std::tuple<StaticArray<Value,Dim+1>>,
        std::tuple<
            StaticArray<Value,Dim+1>,
            TData>
    > _data;
};


template <class T, class TData = void>
concept QuadraturePointLike
= requires (T& instance, const T& constInstance) {
    {constInstance.evaluate(
        maths::IdentityTransform<double,1u>(),
        std::span<double>())};
    {constInstance.position()};
    {constInstance.weight()};
    {constInstance.data()} -> std::same_as<typename VoidSafe<const TData>::Ref>;
    {instance.data()} -> std::same_as<typename VoidSafe<TData>::Ref>;
}; // concept QuadraturePointLike


} // namespace cie::fem
