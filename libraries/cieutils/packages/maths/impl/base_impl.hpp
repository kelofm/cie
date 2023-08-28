#ifndef CIE_UTILS_MATHS_BASE_IMPL_HPP
#define CIE_UTILS_MATHS_BASE_IMPL_HPP

// --- Internal Includes ---
#include "packages/macros/inc/exceptions.hpp"
#include "packages/stl_extension/inc/StaticArray.hpp"

// --- STL Includes ---
#include <string>


namespace cie::utils {


template <  Size NDigits,
            concepts::Integer NT >
constexpr void
baseN(  NT base10,
        NT base,
        StaticArray<NT,NDigits>& output )
{
    CIE_BEGIN_EXCEPTION_TRACING

    for (NT i = 0; i < NDigits; ++i)
    {
        if (base10 == 0)
            output[i]   = 0;
        else
        {
            output[i]   = base10 % base;
            base10     /= base;
        }
    }
    if (base10 != 0)
        throw std::overflow_error( "Error converting to base " + std::to_string(base) );

    CIE_END_EXCEPTION_TRACING
}


template <  Size NDigits,
            concepts::Integer NT >
constexpr StaticArray<NT,NDigits>
baseN(  NT base10,
        NT base )
{
    CIE_BEGIN_EXCEPTION_TRACING

    StaticArray<NT,NDigits> result;
    baseN<NDigits>(base10, base, result);
    return result;

    CIE_END_EXCEPTION_TRACING
}


template <  Size NDigits,
            concepts::Integer NT >
constexpr NT
base10( const StaticArray<NT,NDigits>& base_N,
        NT base )
{
    CIE_BEGIN_EXCEPTION_TRACING

    NT power = 1;
    NT base10 = 0;
    for (NT i = 0; i < NDigits; ++i)
    {
        base10 += power*base_N[i];
        power *= base;
    }
    return base10;

    CIE_END_EXCEPTION_TRACING
}


template <  Size NDigits,
            concepts::Integer NT >
constexpr NT
reinterpretBase(    const StaticArray<NT,NDigits>& base_N,
                    NT newBase )
{
    CIE_BEGIN_EXCEPTION_TRACING

    NT power = 1;
    NT value = 0;
    for (NT i=0; i<NDigits; ++i)
    {
        value += power*base_N[i];
        power *= newBase;
    }
    return value;

    CIE_END_EXCEPTION_TRACING
}


} // namespace cie::utils


#endif