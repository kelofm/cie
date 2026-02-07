#pragma once

// help the language server
#include "packages/numeric/inc/CompositeDomain.hpp"



namespace cie::fem {


template <CellLike TCell>
constexpr CellDomain<TCell>::CellDomain() noexcept
    : _pCell(nullptr)
{}

template <CellLike TCell>
constexpr CellDomain<TCell>::CellDomain(Ref<const TCell> rCell) noexcept
    : _pCell(&rCell)
{}


template <CellLike TCell>
void CellDomain<TCell>::whichSubdomain(
        std::span<const typename TCell::Value> points,
        std::span<DomainData> subdomains) {
    using Value = typename TCell::Value;
    // TODO improve the search to take advantage of multiple points.
    for (std::size_t iSubdomain=0ul; iSubdomain<subdomains.size(); ++iSubdomain) {
        const std::span<const typename TCell::Value,TCell::PhysicalDimension> physicalPoint(
            points.data() + iSubdomain * TCell::PhysicalDimension,
            TCell::PhysicalDimension);
        std::array<Value,TCell::ParametricDimension> parametricPoint;
        _pCell->transform(
            Kernel<TCell::PhysicalDimension,Value>::template cast<PhysicalCoordinate<Value>>(physicalPoint),
            Kernel<TCell::ParametricDimension,Value>::template cast<ParametricCoordinate<Value>>(parametricPoint));
        subdomains[iSubdomain] = std::all_of(
            parametricPoint.begin(),
            parametricPoint.end(),
            [] (Value coordinate) {
                return static_cast<Value>(-1) <= coordinate
                    && coordinate < static_cast<Value>(1);
            });
    } // for iSubdomain in range(subdomains.size())
}


} // namespace cie::fem
