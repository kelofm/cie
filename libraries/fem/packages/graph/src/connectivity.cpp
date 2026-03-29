// --- FEM Includes ---
#include "packages/graph/inc/connectivity.hpp"


namespace cie::fem {


template <unsigned D>
typename AnsatzMap<D>::AnsatzPairs
AnsatzMap<D>::findPairs(
    const OrientedBoundary<D> first,
    const OrientedBoundary<D> second) const noexcept {
        auto it = _topology.find(std::make_pair(first, second));
        AnsatzPairs output {
            .it = it,
            .swap = false};
        if (it != _topology.end() && first != it->first.first)
            output.swap = true;
        return output;
}


template <unsigned D>
std::size_t AnsatzMap<D>::pairCount(AnsatzPairs pairs) const noexcept {
    if (pairs.it != _topology.end()) {
        return pairs.it->second.size();
    } else {
        return 0ul;
    }
}


template <unsigned Dimension>
void AnsatzMap<Dimension>::getPairs(
    AnsatzPairs pairs,
    std::span<std::pair<std::size_t,std::size_t>> output) const {
        if (pairs.it != _topology.end()) {
            CIE_OUT_OF_RANGE_CHECK(pairs.it->second.size() == output.size())
            if (pairs.swap) {
                std::transform(
                    pairs.it->second.begin(),
                    pairs.it->second.end(),
                    output.begin(),
                    [](const auto& rPair) {
                        return std::make_pair(rPair.second, rPair.first);
                    });
            } else {
                std::copy(
                    pairs.it->second.begin(),
                    pairs.it->second.end(),
                    output.begin());
            }
        } else CIE_CHECK(output.size() == 0ul, "")
}


template <unsigned Dimension>
std::size_t AnsatzMap<Dimension>::ansatzCount() const noexcept {
    return _ansatzCount;
}


template <unsigned D>
void makeAnsatzMask(
    std::size_t setSize,
    std::span<std::uint8_t> mask) {
        CIE_CHECK(
            mask.size() == intPow(setSize, D),
            std::format(
                "expecting an ansatz mask with {} components, but got {}",
                intPow(setSize, D), mask.size()))
        std::fill(
            mask.begin(),
            mask.end(),
            static_cast<std::int8_t>(0));
        std::array<std::uint16_t,D> state;
        std::fill(
            state.begin(),
            state.end(),
            static_cast<std::uint8_t>(0));

        std::size_t iAnsatz = 0ul;
        do {
            mask[iAnsatz] = std::max<std::uint8_t>(
                *std::max_element(
                    state.begin(),
                    state.end()),
                1);
            ++iAnsatz;
        } while (cie::maths::OuterProduct<D>::next(setSize, state.data()));
}


#define CIE_INSTANTIATE_ANSATZ_MAP(D)                                       \
    template class AnsatzMap<D>;                                            \
    template void makeAnsatzMask<D>(std::size_t, std::span<std::uint8_t>);


CIE_INSTANTIATE_ANSATZ_MAP(1)
CIE_INSTANTIATE_ANSATZ_MAP(2)
CIE_INSTANTIATE_ANSATZ_MAP(3)


#undef CIE_INSTANTIATE_ANSATZ_MAP


} // namespace cie::fem
