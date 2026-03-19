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


#define CIE_INSTANTIATE_ANSATZ_MAP(D)               \
    template class AnsatzMap<D>;


CIE_INSTANTIATE_ANSATZ_MAP(1)
CIE_INSTANTIATE_ANSATZ_MAP(2)
CIE_INSTANTIATE_ANSATZ_MAP(3)


#undef CIE_INSTANTIATE_ANSATZ_MAP


} // namespace cie::fem
