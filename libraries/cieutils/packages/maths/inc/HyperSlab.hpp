#pragma once

// --- Utility Includes ---
#include "packages/compile_time/packages/concepts/inc/basic_concepts.hpp"

// --- STL Includes ---
#include <array>
#include <span>


namespace cie {


/// @ingroup cieutils
template <unsigned Dim, concepts::Integer TValue = int>
class HyperSlab {
public:
    static inline constexpr unsigned Dimension = Dim;

    constexpr HyperSlab() noexcept;

    explicit constexpr HyperSlab(std::span<const TValue,2*Dim> data) noexcept;

    explicit constexpr HyperSlab(Ref<const std::array<TValue,2*Dim>> rData) noexcept;

    [[nodiscard]] constexpr bool empty() const noexcept;

    [[nodiscard]] constexpr TValue count() const noexcept;

    [[nodiscard]] constexpr HyperSlab intersect(Ref<const HyperSlab> rRhs) const noexcept;

    [[nodiscard]] constexpr bool expand(
        std::span<const TValue,Dim> range,
        std::span<TValue> out) const noexcept;

private:
    std::array<TValue,2*Dimension> _data;
}; // class HyperSlab


} // namespace cie

#include "packages/maths/impl/HyperSlab_impl.hpp"
