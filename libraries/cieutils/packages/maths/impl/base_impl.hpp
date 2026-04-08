#pragma once

// --- Internal Includes ---
#include "packages/stl_extension/inc/StaticArray.hpp"

// --- STL Includes ---
#include <string>


namespace cie::utils {


template <Size NDigits, concepts::Integer NT >
constexpr void baseN(
    NT base10,
    NT base,
    StaticArray<NT,NDigits>& output ) {
        for (NT i = 0; i < NDigits; ++i) {
            if (base10 == 0)
                output[i]   = 0;
            else {
                output[i]   = base10 % base;
                base10     /= base;
            }
        }
        if (base10 != 0)
            throw std::overflow_error( "Error converting to base " + std::to_string(base) );
}


template <Size NDigits, concepts::Integer NT >
constexpr StaticArray<NT,NDigits> baseN(
    NT base10,
    NT base) {
        StaticArray<NT,NDigits> result;
        baseN<NDigits>(base10, base, result);
        return result;
}


template <Size NDigits, concepts::Integer NT >
constexpr NT base10(
    const StaticArray<NT,NDigits>& base_N,
    NT base) {
        NT power = 1;
        NT base10 = 0;
        for (NT i = 0; i < NDigits; ++i) {
            base10 += power*base_N[i];
            power *= base;
        }
        return base10;
}


template <Size NDigits, concepts::Integer NT >
constexpr NT reinterpretBase(
    const StaticArray<NT,NDigits>& base_N,
    NT newBase ) {
        NT power = 1;
        NT value = 0;
        for (NT i=0; i<NDigits; ++i) {
            value += power*base_N[i];
            power *= newBase;
        }
        return value;
}


} // namespace cie::utils
