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





} // namespace cie::fem
