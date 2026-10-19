#pragma once

// --- Utility Includes ---
#include "packages/compile_time/packages/concepts/inc/basic_concepts.hpp"

// --- STL Includes ---
#include <array>
#include <span>
#include <tuple>


namespace cie {


/// @ingroup cieutils
template <unsigned Dim, concepts::Integer TValue = int>
class HyperSlab {
public:
    static inline constexpr unsigned Dimension = Dim;

    using Value = TValue;

    constexpr HyperSlab() noexcept;

    explicit constexpr HyperSlab(std::span<const TValue,2*Dim> data) noexcept;

    explicit constexpr HyperSlab(Ref<const std::array<TValue,2*Dim>> rData) noexcept;

    [[nodiscard]] constexpr std::span<const TValue,2*Dim> data() const noexcept;

    [[nodiscard]] constexpr std::span<TValue,2*Dim> data() noexcept;

    [[nodiscard]] constexpr bool empty() const noexcept;

    [[nodiscard]] constexpr TValue count() const noexcept;

    [[nodiscard]] constexpr HyperSlab intersect(Ref<const HyperSlab> rRhs) const noexcept;

    [[nodiscard]] constexpr bool expand(
        std::span<const TValue,Dim> range,
        std::span<TValue> out) const noexcept;

    [[nodiscard]] constexpr bool expand(
        Ref<const std::array<TValue,Dim>> rRange,
        std::span<TValue> out) const noexcept;

private:
    std::array<TValue,2*Dimension> _data;
}; // class HyperSlab


template <class T>
concept HyperSlabLike = std::is_same_v<
    T,
    HyperSlab<T::Dimension,typename T::Value>>;


template <class T>
concept HyperSlabArrayLike
= requires {
    typename std::remove_cvref_t<T>::value_type;
    typename std::tuple_size<std::remove_cvref_t<T>>::type;
} && HyperSlabLike<typename std::remove_cvref_t<T>::value_type>
&& (std::tuple_size_v<std::remove_cvref_t<T>> > 0ul)
&& std::is_same_v<
    std::remove_cvref_t<T>,
    std::array<
        typename std::remove_cvref_t<T>::value_type,
        std::tuple_size_v<std::remove_cvref_t<T>>>>;


template <class T, unsigned Dim, class TValue>
concept TilingLike
=  std::is_same_v<typename T::Value,HyperSlab<Dim,TValue>>
&& requires (
    T instance,
    HyperSlab<Dim,TValue> target,
    std::span<HyperSlab<Dim,TValue>> outRange) {
        {instance.count(target)}            -> std::same_as<TValue>;
        {instance.tile(target, outRange)}   -> std::same_as<bool>;
}; // concept TilingLike


template <unsigned Dim, class T>
class Tiling {
public:
    using Value = HyperSlab<Dim,T>;

    constexpr Tiling() noexcept;

    constexpr Tiling(Ref<const std::array<T,Dim>> rTileShape) noexcept;

    constexpr Tiling(std::span<const T,Dim> tileShape) noexcept;

    constexpr virtual ~Tiling() = default;

    constexpr virtual T count(Ref<const Value> rTarget) const noexcept;

    constexpr virtual bool tile(
        Ref<const Value> rTarget,
        std::span<Value> out) const noexcept;

private:
    std::array<T,Dim> _tileShape;
}; // class Tiling


} // namespace cie

#include "packages/maths/impl/HyperSlab_impl.hpp"
