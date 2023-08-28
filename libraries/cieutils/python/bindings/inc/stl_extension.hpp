#ifndef CIE_UTILS_STL_EXTENSION_HPP
#define CIE_UTILS_STL_EXTENSION_HPP

// --- External Includes ---
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "pybind11/numpy.h"

// --- Utility Includes ---
#include "packages/stl_extension/inc/StaticArray.hpp"
#include "packages/stl_extension/inc/DynamicArray.hpp"
#include "python/bindings/inc/common.hpp"

// --- STL Includes ---
#include <ranges>
#include <span>


template <class T, class TAllocator>
struct pybind11::detail::type_caster<cie::DynamicArray<T,TAllocator>>
    : pybind11::detail::list_caster<cie::DynamicArray<T,TAllocator>,T>
{};


template <class T, cie::Size ArraySize>
struct pybind11::detail::type_caster<cie::StaticArray<T,ArraySize>>
    : pybind11::detail::array_caster<cie::StaticArray<T,ArraySize>,T,false,ArraySize>
{};


namespace cie {


/// @brief Compute byte strides of a contiguous array.
template <class TValue, unsigned Dimension>
StaticArray<Size,Dimension> makeContiguousStrides(Ptr<const Size> p_shapeBegin)
{
    StaticArray<Size,Dimension> strides;
    const auto shape = std::span(p_shapeBegin, p_shapeBegin + Dimension);
    std::exclusive_scan(shape.rbegin(),
                        shape.rend(),
                        strides.rbegin(),
                        sizeof(TValue),
                        std::multiplies<Size>());
    return strides;
}


/// @brief Wrapper for the @ref pybind11::array_t constructor.
template <class TValue, unsigned Dimension, int Flags = pybind11::array::c_style>
pybind11::array_t<TValue,Flags> makeNumpyArray(TValue* p_begin,
                                               const Size* p_shapeBegin,
                                               const Size* p_strideBegin,
                                               pybind11::capsule sentinel)
{
    return pybind11::array_t<TValue,Flags>(
        std::span(p_shapeBegin, Dimension),
        std::span(p_strideBegin, Dimension),
        p_begin,
        sentinel
    );
}


/// @brief Create a managed contiguous array.
template <class TValue, unsigned Dimension, int Flags = pybind11::array::c_style>
auto makeNumpyArray(Ptr<const Size> p_shapeBegin)
{
    // Allocate a flat array of the exact size
    const Size flatSize = std::accumulate(p_shapeBegin,
                                          p_shapeBegin + Dimension,
                                          static_cast<Size>(1),
                                          std::multiplies<Size>());
    TValue* array = new TValue[flatSize];

    // Define a sentinel that deallocates the array when invoked
    pybind11::capsule release (array, [] (void* p_array) {
        delete[] reinterpret_cast<TValue*>(p_array);
    });

    return makeNumpyArray<TValue,Dimension,Flags>(
        array,
        p_shapeBegin,
        makeContiguousStrides<TValue,Dimension>(p_shapeBegin).data(),
        release
    );
}


} // namespace cie


#endif
