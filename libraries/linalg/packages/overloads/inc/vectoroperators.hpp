#ifndef CIE_LINALG_VECTOR_OPERATORS
#define CIE_LINALG_VECTOR_OPERATORS

// --- Utility Includes ---
#include "packages/compile_time/packages/concepts/inc/container_concepts.hpp"

// --- Internal Includes ---
#include "packages/types/inc/vectortypes.hpp"
#include "packages/types/inc/arraytypes.hpp"


// ---------------------------------------------------------
// VECTOR - SCALAR OPERATORS
// ---------------------------------------------------------

template <cie::concepts::NumericContainer ArrayType, cie::concepts::Numeric ScalarType>
void operator*=( ArrayType& vector, const ScalarType& scalar );

template <cie::concepts::NumericContainer ArrayType, cie::concepts::Numeric ScalarType>
void operator/=( ArrayType& vector, const ScalarType& scalar );

template <cie::concepts::NumericContainer ArrayType, cie::concepts::Numeric ScalarType>
ArrayType operator+( const ArrayType& vector, const ScalarType& scalar );

template <cie::concepts::NumericContainer ArrayType, cie::concepts::Numeric ScalarType>
ArrayType operator+( const ScalarType& scalar, const ArrayType& vector );

template <cie::concepts::NumericContainer ArrayType, cie::concepts::Numeric ScalarType>
ArrayType operator-( const ArrayType& vector, const ScalarType& scalar );

template <cie::concepts::NumericContainer ArrayType, cie::concepts::Numeric ScalarType>
ArrayType operator-( const ScalarType& scalar, const ArrayType& vector );

template <cie::concepts::NumericContainer ArrayType, cie::concepts::Numeric ScalarType>
ArrayType operator*( const ArrayType& vector, const ScalarType& scalar );

template <cie::concepts::NumericContainer ArrayType, cie::concepts::Numeric ScalarType>
ArrayType operator*( const ScalarType& scalar, const ArrayType& vector );

template <cie::concepts::NumericContainer ArrayType, cie::concepts::Numeric ScalarType>
ArrayType operator/( const ArrayType& vector, const ScalarType& scalar );


// ---------------------------------------------------------
// VECTOR - VECTOR OPERATORS
// ---------------------------------------------------------

template <cie::concepts::NumericContainer ArrayType>
void operator+=( ArrayType& r_lhs, const ArrayType& r_rhs );

template <cie::concepts::NumericContainer ArrayType>
void operator-=( ArrayType& r_lhs, const ArrayType& r_rhs );

template <cie::concepts::NumericContainer ArrayType>
ArrayType operator+( const ArrayType& lhs, const ArrayType& rhs );

template <cie::concepts::NumericContainer ArrayType>
ArrayType operator-( const ArrayType& lhs, const ArrayType& rhs );

template <cie::concepts::NumericContainer ArrayType>
typename ArrayType::value_type operator*( const ArrayType& lhs, const ArrayType& rhs );


#include "packages/overloads/impl/vectoroperators_impl.hpp"

#endif