#ifndef LINALG_TYPE_OPERATIONS_HPP
#define LINALG_TYPE_OPERATIONS_HPP

// --- Utility Includes ---
#include "packages/compile_time/packages/concepts/inc/container_concepts.hpp"

// --- Internal Includes ---
#include "packages/types/inc/functiontypes.hpp"

namespace cie::linalg {

template <concepts::NumericContainer T>
typename T::value_type norm2( const T& r_object );

template <concepts::NumericContainer T>
typename T::value_type norm( const T& r_object );

template <concepts::NumericContainer T>
void normalize( T& r_object,
                NormFunction<T> measure=&norm<T>,
                typename T::value_type tolerance=1e-15 );

} // namespace cie::linalg

#include "packages/types/impl/typeoperations_impl.hpp"

#endif