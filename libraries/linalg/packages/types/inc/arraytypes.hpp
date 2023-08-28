#ifndef LINALG_ARRAY_TYPES_HPP
#define LINALG_ARRAY_TYPES_HPP

// --- Utility Includes ---
#include "packages/stl_extension/inc/StaticArray.hpp"

// --- STL Includes ---
#include <array>
#include <cstddef>
#include <cstdint>

namespace cie {

typedef StaticArray<double, 2>   Point2;
typedef StaticArray<double, 3>   Point3;

template <size_t N>
using DoubleArray = StaticArray<double, N>;

template <size_t N>
using FloatArray = StaticArray<float, N>;

template <size_t N>
using IntArray = StaticArray<int, N>;

template <size_t N>
using UIntArray     = StaticArray<size_t, N>;
template <size_t N>
using UInt8Array    = StaticArray<uint8_t, N>;

template <size_t N>
using UInt16Array   = StaticArray<uint16_t, N>;

} // namespace cie

#endif