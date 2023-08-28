#ifndef CIE_UTILS_MATHS_N_CHOOSE_K_HPP
#define CIE_UTILS_MATHS_N_CHOOSE_K_HPP

// --- Utility Inlcudes ---
#include "packages/types/inc/types.hpp"

// --- STL Includes ---
#include <vector>
#include <string>


namespace cie::utils {


class NChooseK
{
public:
    using state_container = std::vector<Size>;

public:
    NChooseK( Size n, Size k );

    NChooseK() = delete;

    NChooseK( NChooseK&& r_rhs ) = delete;

    NChooseK( const NChooseK& r_rhs ) = delete;

    NChooseK& operator=( const NChooseK& r_rhs ) = delete;

    bool next();

    void reset();

    state_container::const_iterator begin() const;

    state_container::const_iterator end() const;

    Size n() const;

    Size k() const;

    Size numberOfPermutations() const;

private:
    static Size nChooseK( Size n, Size k );

private:
    state_container _state;
    std::string     _mask;
};


} // namespace cie::utils


#endif