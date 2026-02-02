#pragma once

// --- FEM Includes ---
#include "packages/numeric/inc/Cell.hpp"

// --- STL Includes ---
#include <span>
#include <type_traits>


namespace cie::fem {


template <class T>
concept DomainDataLike
=  std::is_default_constructible_v<T>
&& requires (T instance) {
    {instance == instance}  -> std::same_as<bool>;
    {instance != instance}  -> std::same_as<bool>;
    {bool(instance)}        -> std::same_as<bool>;
}; // concept DomainDataLike


template <class T, unsigned Dimension, class TValue>
concept CompositeDomainLike
=  DomainDataLike<typename T::DomainData>
&& requires (const T& constInstance,
             std::span<const TValue> points,
             std::span<typename T::DomainData> subdomains) {
    {constInstance.whichSubdomain(points, subdomains)};
}; // CompositeDomain


template <CellLike TCell>
class CellDomain {
public:
    using Cell = TCell;

    using DomainData = bool;

    constexpr CellDomain() noexcept;

    constexpr CellDomain(Ref<const TCell> rCell) noexcept;

    void whichSubdomain(
        std::span<const typename TCell::Value> points,
        std::span<DomainData> subdomains);

private:
    Ptr<const TCell> _pCell;
}; // class CellDomain


} // namespace cie::fem

#include "packages/numeric/impl/CompositeDomain_impl.hpp"
