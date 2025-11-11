// --- Utility Includes ---
#include "packages/macros/inc/exceptions.hpp"
#include "packages/macros/inc/checks.hpp"
#include "packages/stl_extension/inc/Hash.hpp"
#include "packages/maths/inc/NChooseK.hpp"

// --- STL Includes ---
#include <numeric>
#include <algorithm>
#include <unordered_map>


namespace cie::utils {


Size nChooseK(unsigned n,
              unsigned k,
              Ref<std::unordered_map<
                std::pair<unsigned,unsigned>,
                Size,
                Hash<std::pair<unsigned,unsigned>>,
                Equal<std::pair<unsigned,unsigned>>
              >> rMemo) noexcept
{
    // Special cases.
    if      (k == 0u || n == k) return 1ul;
    else if (k == 1u || k == n - 1u) return static_cast<Size>(n);

    // From this point 0 < k && 0 < n.

    // Check for an entry already computed.
    auto emplaceResult = rMemo.emplace(std::pair<unsigned,unsigned> {n, k},
                            0ul);

    if (emplaceResult.second) {
        // No entry found => do recursion.
        emplaceResult.first->second = k + k < n
                                      ? n * nChooseK(n - 1u, k - 1u, rMemo) / k
                                      : n * nChooseK(n - 1u, k,      rMemo) / (n - k);
    }

    return emplaceResult.first->second;
}


struct NChooseK::Impl {
    state_container _state;

    std::string     _mask;

    std::unordered_map<
        std::pair<unsigned,unsigned>,
        Size,
        Hash<std::pair<unsigned,unsigned>>,
        Equal<std::pair<unsigned,unsigned>>
    > _memo;
}; // struct NChooseK::Impl


NChooseK::NChooseK(unsigned n, unsigned k)
    : _pImpl(new Impl {
        ._state = state_container(k),
        ._mask = std::string(n, 0),
        ._memo = {}
    })
{
    CIE_CHECK(k <= n, "Cannot choose " + std::to_string(k) + " from " + std::to_string(n) + " elements!")

    CIE_BEGIN_EXCEPTION_TRACING
    this->reset();
    CIE_END_EXCEPTION_TRACING
}


NChooseK::~NChooseK() = default;


bool NChooseK::next() noexcept
{
    auto& rMask = this->_pImpl->_mask;
    auto& rState = this->_pImpl->_state;

    if (!std::prev_permutation(rMask.begin(), rMask.end()))
        return false;

    auto it_state = rState.begin();

    for (Size i=0; i<this->n(); ++i)
        if (rMask[i])
            *it_state++ = i;

    return true;
}


void NChooseK::reset() noexcept
{
    std::iota(
        this->_pImpl->_state.begin(),
        this->_pImpl->_state.end(),
        0
    );

    std::fill_n(
        this->_pImpl->_mask.begin(),
        this->k(),
        1
    );

    std::fill(
        this->_pImpl->_mask.begin() + this->k(),
        this->_pImpl->_mask.end(),
        0
    );
}


NChooseK::state_container::const_iterator NChooseK::begin() const noexcept
{
    return this->_pImpl->_state.begin();
}


NChooseK::state_container::const_iterator NChooseK::end() const noexcept
{
    return this->_pImpl->_state.end();
}


unsigned NChooseK::n() const noexcept
{
    return static_cast<unsigned>(this->_pImpl->_mask.size());
}


unsigned NChooseK::k() const noexcept
{
    return static_cast<unsigned>(this->_pImpl->_state.size());
}


Size NChooseK::numberOfPermutations() const noexcept
{
    return nChooseK(this->n(), this->k(), this->_pImpl->_memo);
}


} // namespace cie::utils
