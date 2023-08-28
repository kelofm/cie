#ifndef LINALG_TYPE_OPERATIONS_IMPL_HPP
#define LINALG_TYPE_OPERATIONS_IMPL_HPP

// --- Utility Includes ---
#include "packages/macros/inc/exceptions.hpp"

// --- Linalg Includes ---
#include "packages/types/inc/functiontypes.hpp"
#include "packages/types/inc/typeoperations.hpp"

// --- STL Includes ---
#include <numeric>
#include <cmath>
#include <stdexcept>


namespace cie::linalg {

template <concepts::NumericContainer T>
inline typename T::value_type
norm2( const T& r_object )
{
    return std::accumulate(
        r_object.begin(),
        r_object.end(),
        0.0,
        []( auto lhs, auto rhs ) { return lhs + rhs*rhs; }
    );
}


template <concepts::NumericContainer T>
inline typename T::value_type
norm( const T& r_object )
{
    return std::sqrt( norm2<T>(r_object) );
}


template <concepts::NumericContainer T>
inline void
normalize(T& r_object,
          NormFunction<T> measure,
          typename T::value_type tolerance)
{
    typename T::value_type objectNorm = measure(r_object);

    if ( std::abs(objectNorm)<tolerance )
        CIE_THROW( std::runtime_error, "Cannot normalize an object with 0 norm!" )
    else
        for ( auto& r_component : r_object )
            r_component /= objectNorm;
}

} // namespace cie::linalg

#endif