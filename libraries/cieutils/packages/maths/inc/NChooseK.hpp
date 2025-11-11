#pragma once

// --- Utility Inlcudes ---
#include "packages/types/inc/types.hpp"

// --- STL Includes ---
#include <vector>
#include <string>
#include <memory>


namespace cie::utils {


/// @ingroup cieutils
class NChooseK
{
public:
    using state_container = std::vector<unsigned>;

public:
    NChooseK(unsigned n, unsigned k);

    NChooseK() = delete;

    NChooseK(NChooseK&& r_rhs) = delete;

    NChooseK(const NChooseK& r_rhs) = delete;

    NChooseK& operator=(const NChooseK& r_rhs) = delete;

    ~NChooseK();

    bool next() noexcept;

    void reset() noexcept;

    state_container::const_iterator begin() const noexcept;

    state_container::const_iterator end() const noexcept;

    unsigned n() const noexcept;

    unsigned k() const noexcept;

    Size numberOfPermutations() const noexcept;

private:
    struct Impl;
    std::unique_ptr<Impl> _pImpl;
};


} // namespace cie::utils
