#pragma once

// --- FEM Includes ---
#include "packages/numeric/inc/MeshBase.hpp"

// --- STL Includes ---
#include <algorithm>


namespace cie::fem {


template <maths::Expression TAnsatzSpace>
MeshBase<TAnsatzSpace>::MeshBase(std::span<const TAnsatzSpace> spaces) {
    CIE_BEGIN_EXCEPTION_TRACING
        _spaces.reserve(spaces.size());
        std::transform(
            spaces.begin(),
            spaces.end(),
            std::back_inserter(_spaces),
            [] (Ref<const TAnsatzSpace> rSpace) -> typename decltype(_spaces)::value_type {
                return {
                    rSpace,
                    rSpace.makeDerivative()};
            });
    CIE_END_EXCEPTION_TRACING
}


template <maths::Expression TAnsatzSpace>
constexpr std::size_t MeshBase<TAnsatzSpace>::ansatzCount() const noexcept {
    return _spaces.size();
}


template <maths::Expression TAnsatzSpace>
Ref<const TAnsatzSpace>
MeshBase<TAnsatzSpace>::ansatz(std::size_t iAnsatz) const noexcept {
    return std::get<0>(_spaces[iAnsatz]);
}


template <maths::Expression TAnsatzSpace>
Ref<const typename TAnsatzSpace::Derivative>
MeshBase<TAnsatzSpace>::ansatzDerivative(std::size_t iAnsatz) const noexcept {
    return std::get<1>(_spaces[iAnsatz]);
}

template <maths::Expression TAS>
void io::GraphML::Deserializer<MeshBase<TAS>>::onElementBegin(
    [[maybe_unused]] Ptr<void> pThis,
    [[maybe_unused]] std::string_view elementName,
    [[maybe_unused]] std::span<io::GraphML::AttributePair> attributes) {
        CIE_THROW(NotImplementedException, "")
}

template <maths::Expression TAS>
void io::GraphML::Deserializer<MeshBase<TAS>>::onText(
    [[maybe_unused]] Ptr<void> pThis,
    [[maybe_unused]] std::string_view elementName) {
        CIE_THROW(NotImplementedException, "")
}

template <maths::Expression TAS>
void io::GraphML::Deserializer<MeshBase<TAS>>::onElementEnd(
    [[maybe_unused]] Ptr<const void> pThis,
    [[maybe_unused]] std::string_view elementName) {
        CIE_THROW(NotImplementedException, "")
}


} // namespace cie::fem
