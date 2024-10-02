#ifndef CIE_FEM_ORIENTED_BOUNDARY_IMPL_HPP
#define CIE_FEM_ORIENTED_BOUNDARY_IMPL_HPP

// --- Utility Includes ---
#include "packages/maths/inc/power.hpp"
#include "packages/macros/inc/checks.hpp"
#include "packages/macros/inc/exceptions.hpp"

// Help the language server
#include "packages/graph/inc/OrientedBoundary.hpp"


namespace cie::fem {


template <unsigned Dimension>
template <bool Mutable>
OrientedBoundary<Dimension>::ValueProxy<Mutable>::ValueProxy(Ref<std::conditional_t<Mutable,Data,const Data>> rData,
                                                             unsigned index)
    : _pData(&rData),
      _index(index)
{
}


template <unsigned Dimension>
template <bool Mutable>
typename OrientedBoundary<Dimension>::template ValueProxy<Mutable>&
OrientedBoundary<Dimension>::ValueProxy<Mutable>::operator=(BoundaryID rhs)
requires Mutable
{
    Ref<Data> rData = *_pData;

    // Zero out the region to be overwritten.
    const unsigned bitOffset = _index * ComponentBitWidth;
    constexpr Data localMask = intPow(2u, ComponentBitWidth) - 1;
    const Data mask = localMask << bitOffset;
    rData &= ~mask;

    // Compute the new component.
    Data componentData = rhs.getDimension();
    componentData <<= 1u;
    componentData |= Data(rhs.getDirection());
    componentData <<= bitOffset;

    // Overwrite the component.
    rData |= componentData;

    return *this;
}


template <unsigned Dimension>
template <bool Mutable>
OrientedBoundary<Dimension>::ValueProxy<Mutable>::operator BoundaryID() const noexcept
{
    constexpr Data mask = intPow(2u, ComponentBitWidth) - 1u;
    const unsigned bitOffset = _index * ComponentBitWidth;

    // Get bitwise data to the least significant region of an integer.
    Ref<const Data> rData = *_pData;
    const unsigned componentData = ((rData >> bitOffset) & mask).to_ulong();

    // The first bit indicates positive/negative direction, while the rest
    // define an integer that encodes the index of the axis.
    const bool direction = componentData & 1u;
    const unsigned axis = componentData >> 1;

    return BoundaryID(axis, direction);
}


template <unsigned Dimension>
template <bool Mutable>
OrientedBoundary<Dimension>::Iterator<Mutable>::Iterator(Iterator<true>&& rhs) noexcept
requires (!Mutable)
    : Iterator(rhs._pData, rhs._index)
{
}


template <unsigned Dimension>
template <bool Mutable>
OrientedBoundary<Dimension>::Iterator<Mutable>::Iterator(const Iterator<true>& rhs) noexcept
requires (!Mutable)
    : Iterator(rhs._pData, rhs._index)
{
}


template <unsigned Dimension>
template <bool Mutable>
OrientedBoundary<Dimension>::Iterator<Mutable>::Iterator(Ref<std::conditional_t<Mutable,Data,const Data>> rData,
                                                         unsigned index) noexcept
    : _pData(&rData),
      _index(index)
{
}


template <unsigned Dimension>
template <bool Mutable>
typename OrientedBoundary<Dimension>::template Iterator<Mutable>::reference
OrientedBoundary<Dimension>::Iterator<Mutable>::operator*() const
{
    return typename OrientedBoundary<Dimension>::template ValueProxy<Mutable>(*_pData, _index);
}


template <unsigned Dimension>
template <bool Mutable>
typename OrientedBoundary<Dimension>::template Iterator<Mutable>&
OrientedBoundary<Dimension>::Iterator<Mutable>::operator++() noexcept
{
    ++_index;
    return *this;
}


template <unsigned Dimension>
template <bool Mutable>
typename OrientedBoundary<Dimension>::template Iterator<Mutable>
OrientedBoundary<Dimension>::Iterator<Mutable>::operator++(int) noexcept
{
    Iterator copy = *this;
    ++(*this);
    return copy;
}


template <unsigned Dimension>
template <bool Mutable>
typename OrientedBoundary<Dimension>::template Iterator<Mutable>&
OrientedBoundary<Dimension>::Iterator<Mutable>::operator+=(difference_type offset) noexcept
{
    _index += offset;
    return *this;
}


template <unsigned Dimension>
template <bool Mutable>
typename OrientedBoundary<Dimension>::template Iterator<Mutable>
OrientedBoundary<Dimension>::Iterator<Mutable>::operator+(difference_type offset) noexcept
{
    Iterator copy = *this;
    copy += offset;
    return copy;
}


template <unsigned Dimension>
template <bool Mutable>
typename OrientedBoundary<Dimension>::template Iterator<Mutable>&
OrientedBoundary<Dimension>::Iterator<Mutable>::operator--() noexcept
{
    --_index;
    return *this;
}


template <unsigned Dimension>
template <bool Mutable>
typename OrientedBoundary<Dimension>::template Iterator<Mutable>
OrientedBoundary<Dimension>::Iterator<Mutable>::operator--(int) noexcept
{
    Iterator copy = *this;
    --(*this);
    return copy;
}

template <unsigned Dimension>
template <bool Mutable>
typename OrientedBoundary<Dimension>::template Iterator<Mutable>&
OrientedBoundary<Dimension>::Iterator<Mutable>::operator-=(difference_type offset) noexcept
{
    _index -= offset;
    return *this;
}


template <unsigned Dimension>
template <bool Mutable>
typename OrientedBoundary<Dimension>::template Iterator<Mutable>
OrientedBoundary<Dimension>::Iterator<Mutable>::operator-(difference_type offset) noexcept
{
    Iterator copy = *this;
    copy -= offset;
    return copy;
}


template <unsigned Dimension>
template <bool Mutable>
typename OrientedBoundary<Dimension>::template Iterator<Mutable>::difference_type
OrientedBoundary<Dimension>::Iterator<Mutable>::operator-(Iterator rhs) noexcept
{
    return _index - rhs._index;
}


template <unsigned Dimension>
template <bool Mutable>
template <bool M>
bool OrientedBoundary<Dimension>::Iterator<Mutable>::operator==(Iterator<M> right) noexcept
{
    return _pData == right._pData && _index == right._index;
}


template <unsigned Dimension>
template <bool Mutable>
template <bool M>
bool OrientedBoundary<Dimension>::Iterator<Mutable>::operator!=(Iterator<M> right) noexcept
{
    return _index != right._index || _pData != right._pData;
}


template <unsigned Dimension>
OrientedBoundary<Dimension>::OrientedBoundary() noexcept
    : _data(0u)
{
    // By default, the system is not rotated, so all axes align with
    // the original ones.
    for (unsigned iDimension=0u; iDimension<Dimension; ++iDimension) {
        this->at(iDimension) = BoundaryID(iDimension, true);
    }

    // Default boundary refers to the negative side of
    // the first direction ("-x").
    this->id() = BoundaryID();
}


template <unsigned Dimension>
OrientedBoundary<Dimension>::OrientedBoundary(Ptr<const BoundaryID> itBegin,
                                              size_type size,
                                              BoundaryID id)
{
    CIE_OUT_OF_RANGE_CHECK(size == this->size())
    std::copy(itBegin, itBegin + size, this->begin());
    this->id() = id;
}


template <unsigned Dimension>
OrientedBoundary<Dimension>::OrientedBoundary(OrientedAxes<Dimension> axes,
                                              BoundaryID id)
{
    auto it = this->begin();
    for (const auto axis : axes) *it++ = axis;
    this->id() = id;
}


template <unsigned Dimension>
OrientedBoundary<Dimension>::OrientedBoundary(const char axes[2 * Dimension + 1],
                                              const char id[3])
requires (Dimension < 4)
{
    char axis[3];
    axis[2] = '\0';

    for (unsigned iChar=0u; iChar<2*Dimension; iChar+=2) {
        axis[0] = axes[iChar];
        axis[1] = axes[iChar + 1];

        CIE_BEGIN_EXCEPTION_TRACING
        this->at(iChar / 2) = BoundaryID(axis);
        CIE_END_EXCEPTION_TRACING
    }

    CIE_BEGIN_EXCEPTION_TRACING
    this->id() = id;
    CIE_END_EXCEPTION_TRACING
}


template <unsigned Dimension>
BoundaryID OrientedBoundary<Dimension>::id() const noexcept
{
    // The boundary ID is stored right after the array of axis IDs,
    // so this is not a segfault.
    return *this->end();
}


template <unsigned Dimension>
typename OrientedBoundary<Dimension>::reference
OrientedBoundary<Dimension>::id() noexcept
{
    return *this->end();
}


template <unsigned Dimension>
BoundaryID OrientedBoundary<Dimension>::localID() const noexcept
{
    const BoundaryID id = this->id();
    const BoundaryID axis = this->at(id.getDimension());
    const bool direction = id.getDirection() == axis.getDirection();
    return BoundaryID(axis.getDimension(), direction);
}


template <unsigned Dimension>
bool OrientedBoundary<Dimension>::operator==(OrientedBoundary rhs) const noexcept
{
    return this->_data == rhs._data;
}


template <unsigned Dimension>
bool OrientedBoundary<Dimension>::operator!=(OrientedBoundary rhs) const noexcept
{
    return this->_data != rhs._data;
}


template <unsigned Dimension>
bool OrientedBoundary<Dimension>::operator<(OrientedBoundary rhs) const noexcept
{
    {
        const auto left = this->id();
        const auto right = rhs.id();

        if (left < right) {
            return true;
        } else if (right < left) {
            return false;
        }
    }

    for (unsigned iAxis=0u; iAxis<this->size(); ++iAxis) {
        const auto left = this->at(iAxis);
        const auto right = rhs[iAxis];

        if (left < right) {
            return true;
        } else if (right < left) {
            return false;
        }
    } // for iAxis in range(this->size)

    return false;
}


template <unsigned Dimension>
OrientedBoundary<Dimension>&
OrientedBoundary<Dimension>::operator++() noexcept
{
    ++this->id();
    return *this;
}


template <unsigned Dimension>
OrientedBoundary<Dimension>
OrientedBoundary<Dimension>::operator++(int) noexcept
{
    const auto copy = *this;
    ++this->id();
    return copy;
}


template <unsigned Dimension>
OrientedBoundary<Dimension>
OrientedBoundary<Dimension>::operator-() const noexcept
{
    OrientedBoundary output = *this;
    output.id() = -BoundaryID(output.id());
    return output;
}


template <unsigned Dimension>
typename OrientedBoundary<Dimension>::const_reference
OrientedBoundary<Dimension>::operator[](size_type index) const
{
    return *(this->begin() + index);
}


template <unsigned Dimension>
typename OrientedBoundary<Dimension>::reference
OrientedBoundary<Dimension>::operator[](size_type index)
{
    return *(this->begin() + index);
}


template <unsigned Dimension>
typename OrientedBoundary<Dimension>::const_reference
OrientedBoundary<Dimension>::at(size_type index) const
{
    return *(this->begin() + index);
}


template <unsigned Dimension>
typename OrientedBoundary<Dimension>::reference
OrientedBoundary<Dimension>::at(size_type index)
{
    return *(this->begin() + index);
}


template <unsigned Dimension>
typename OrientedBoundary<Dimension>::const_iterator
OrientedBoundary<Dimension>::begin() const noexcept
{
    return const_iterator(_data, 0u);
}


template <unsigned Dimension>
typename OrientedBoundary<Dimension>::iterator
OrientedBoundary<Dimension>::begin() noexcept
{
    return iterator(_data, 0u);
}


template <unsigned Dimension>
typename OrientedBoundary<Dimension>::const_iterator
OrientedBoundary<Dimension>::end() const noexcept
{
    return const_iterator(_data, Dimension);
}


template <unsigned Dimension>
typename OrientedBoundary<Dimension>::iterator
OrientedBoundary<Dimension>::end() noexcept
{
    return iterator(_data, Dimension);
}


} // namespace cie::fem


#endif
