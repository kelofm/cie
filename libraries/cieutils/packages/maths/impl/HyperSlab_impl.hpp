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
    :   HyperSlab(std::span<const T,2*D>(
            rData.data(),
            rData.size()))
{}


template <unsigned D, concepts::Integer T>
constexpr std::span<const T,2*D> HyperSlab<D,T>::data() const noexcept {
    return _data;
}


template <unsigned D, concepts::Integer T>
constexpr std::span<T,2*D> HyperSlab<D,T>::data() noexcept {
    return _data;
}


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


template <unsigned D, concepts::Integer T>
constexpr bool HyperSlab<D,T>::expand(
    Ref<const std::array<T,D>> rRange,
    std::span<T> out) const noexcept {
        return this->expand(
            std::span<const T,D>(rRange),
            out);
}


template <unsigned D, class T>
constexpr Tiling<D,T>::Tiling() noexcept
    : _tileShape() {
        std::fill(
            _tileShape.begin(),
            _tileShape.end(),
            static_cast<T>(1));
}


template <unsigned D, class T>
constexpr Tiling<D,T>::Tiling(Ref<const std::array<T,D>> rTileShape) noexcept
    : Tiling(std::span<const T,D>(rTileShape.data(), rTileShape.size()))
{}


template <unsigned D, class T>
constexpr Tiling<D,T>::Tiling(std::span<const T,D> tileShape) noexcept
    : _tileShape() {
        std::transform(
            tileShape.begin(),
            tileShape.end(),
            _tileShape.begin(),
            [] (T component) {
                return std::max<T>(component, 1);
            });
}


template <unsigned D, class T>
constexpr T Tiling<D,T>::count(Ref<const Value> rTarget) const noexcept {
    T out = 1;
    for (unsigned iDimension=0u; iDimension<D; ++iDimension)
        out *= (rTarget.data()[2 * iDimension + 1] - rTarget.data()[2 * iDimension] + _tileShape[iDimension] - 1) / _tileShape[iDimension];
    return out;
}


template <unsigned D, class T>
constexpr bool Tiling<D,T>::tile(
    Ref<const Value> rTarget,
    std::span<Value> out) const noexcept {
        if (out.size() != this->count(rTarget))
            return false;

        std::array<T,D> tileCounts;
        for (unsigned iDimension=0u; iDimension<D; ++iDimension) {
            const T extent = rTarget.data()[2 * iDimension + 1] - rTarget.data()[2 * iDimension];
            tileCounts[iDimension] = (extent + _tileShape[iDimension] - 1) / _tileShape[iDimension];
        }

        std::array<T,D> tileIndices;
        std::fill(
            tileIndices.begin(),
            tileIndices.end(),
            static_cast<T>(0));

        std::array<T,2 * D> tileData;
        for (std::size_t iOut=0; iOut<out.size(); ++iOut) {
            for (unsigned iDimension=0u; iDimension<D; ++iDimension) {
                const T targetBegin = rTarget.data()[2 * iDimension];
                const T targetEnd = rTarget.data()[2 * iDimension + 1];
                const T tileBegin = targetBegin + tileIndices[iDimension] * _tileShape[iDimension];
                const T tileExtent = std::min<T>(
                    _tileShape[iDimension],
                    targetEnd - tileBegin);
                tileData[2 * iDimension] = tileBegin;
                tileData[2 * iDimension + 1] = tileBegin + tileExtent;
            }
            out[iOut] = Value(tileData);

            for (unsigned iDimension=0u; iDimension<D; ++iDimension) {
                if (++tileIndices[iDimension] < tileCounts[iDimension])
                    break;
                tileIndices[iDimension] = 0;
            }
        } // for iOut
        return true;
}


} // namespace cie::fem
