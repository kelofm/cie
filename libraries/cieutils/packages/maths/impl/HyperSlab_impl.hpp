#pragma once

// --- Utility Includes ---
#include "packages/maths/inc/HyperSlab.hpp"

// --- STL Includes ---
#include <algorithm>
#include <type_traits>


namespace cie {


template <unsigned D, concepts::Integer T>
constexpr HyperSlab<D,T>::HyperSlab() noexcept {
    std::fill(
        _data.begin(),
        _data.end(),
        static_cast<T>(0));
}


template <unsigned D, concepts::Integer T>
constexpr HyperSlab<D,T>::HyperSlab(std::span<const T,2*D> data) noexcept {
    std::copy(
        data.begin(),
        data.end(),
        _data.begin());
    for (auto it=_data.data(); it!=_data.end(); it+=2)
        if (*(it + 1) < *it)
            std::swap(*it, *(it + 1));
}


template <unsigned D, concepts::Integer T>
constexpr HyperSlab<D,T>::HyperSlab(Ref<const std::array<T,2*D>> rData) noexcept
    : HyperSlab(std::span<const T,2*D>(
        rData.data(),
        2 * rData.size()))
{}


template <unsigned D, concepts::Integer T>
constexpr bool HyperSlab<D,T>::empty() const noexcept {
    for (auto it=_data.data(); it!=_data.end(); it+=2)
        if (*it == *(it + 1)) return true;
    return false;
}


template <unsigned D, concepts::Integer T>
constexpr T HyperSlab<D,T>::count() const noexcept {
    T out = 1;
    for (auto it=_data.data(); it!=_data.end(); it+=2)
        out *= *(it + 1) - *it;
    return out;
}


template <unsigned D, concepts::Integer T>
constexpr HyperSlab<D,T> HyperSlab<D,T>::intersect(Ref<const HyperSlab<D,T>> rRhs) const noexcept {
    HyperSlab<D,T> out;
    Ptr<const T> it_lhs = _data.data(), it_rhs = rRhs._data.data();
    Ptr<T> it_out = out._data().data();
    for (; it_lhs!=_data.data(); it_lhs+=2, it_rhs+=2, it_out+=2) {
        *it_out = std::max<T>(
            *it_lhs,
            *it_rhs);
        *(it_out + 1) = std::min<T>(
            *(it_lhs + 1),
            *(it_rhs + 1));
        if (*(it_out + 1) < *it_out)
            std::swap(
                *it_out,
                *(it_out + 1));
    }
    return out;
}


template <unsigned D, concepts::Integer T>
constexpr bool HyperSlab<D,T>::expand(
    std::span<const T,D> range,
    std::span<T> out) const noexcept {
        const T entryCount = this->count();
        if (out.size() < static_cast<std::size_t>(entryCount)) return false;

        std::array<T,D> index;
        std::array<T,D> stride;
        for (unsigned dimension=0; dimension<D; ++dimension) {
            if (_data[2*dimension + 1] > range[dimension])
                    return false;
            if constexpr (std::is_signed_v<T>) {
                if (range[dimension] < 0 || _data[2*dimension] < 0)
                    return false;
            }
            index[dimension] = _data[2*dimension];
        }

        T currentStride = 1;
        for (unsigned dimension=0; dimension<D; ++dimension) {
            stride[dimension] = currentStride;
            currentStride *= range[dimension];
        }

        for (T iOut=0; iOut<entryCount; ++iOut) {
            T iFlat = 0;
            for (unsigned dimension=0; dimension<D; ++dimension)
                iFlat += index[dimension] * stride[dimension];
            out[iOut] = iFlat;

            for (unsigned dimension=0; dimension<D; ++dimension) {
                if (++index[dimension] < _data[2*dimension + 1]) break;
                index[dimension] = _data[2*dimension];
            }
        }
        return true;
}


} // namespace cie::fem
