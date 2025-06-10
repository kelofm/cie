#ifndef CIE_FEM_ORIENTED_AXES_IMPL_HPP
#define CIE_FEM_ORIENTED_AXES_IMPL_HPP

// --- Utility Includes ---
#include "packages/maths/inc/power.hpp"
#include "packages/macros/inc/checks.hpp"
#include "packages/macros/inc/exceptions.hpp"

// Help the language server
#include "packages/graph/inc/OrientedAxes.hpp"


namespace cie::fem {


template <unsigned Dimension>
template <bool Mutable>
OrientedAxes<Dimension>::ValueProxy<Mutable>::ValueProxy(Ref<std::conditional_t<Mutable,Data,const Data>> rData,
                                                         unsigned index)
    : _pData(&rData),
      _index(index)
{
}


template <unsigned Dimension>
template <bool Mutable>
typename OrientedAxes<Dimension>::template ValueProxy<Mutable>&
OrientedAxes<Dimension>::ValueProxy<Mutable>::operator=(BoundaryID rhs)
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
OrientedAxes<Dimension>::ValueProxy<Mutable>::operator BoundaryID() const noexcept
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
OrientedAxes<Dimension>::Iterator<Mutable>::Iterator(Iterator<true>&& rhs) noexcept
requires (!Mutable)
    : Iterator(rhs._pData, rhs._index)
{
}


template <unsigned Dimension>
template <bool Mutable>
OrientedAxes<Dimension>::Iterator<Mutable>::Iterator(const Iterator<true>& rhs) noexcept
requires (!Mutable)
    : Iterator(rhs._pData, rhs._index)
{
}


template <unsigned Dimension>
template <bool Mutable>
OrientedAxes<Dimension>::Iterator<Mutable>::Iterator(Ref<std::conditional_t<Mutable,Data,const Data>> rData,
                                                     unsigned index) noexcept
    : _pData(&rData),
      _index(index)
{
}


template <unsigned Dimension>
template <bool Mutable>
typename OrientedAxes<Dimension>::template Iterator<Mutable>::reference
OrientedAxes<Dimension>::Iterator<Mutable>::operator*() const
{
    return typename OrientedAxes<Dimension>::template ValueProxy<Mutable>(*_pData, _index);
}


template <unsigned Dimension>
template <bool Mutable>
typename OrientedAxes<Dimension>::template Iterator<Mutable>&
OrientedAxes<Dimension>::Iterator<Mutable>::operator++() noexcept
{
    ++_index;
    return *this;
}


template <unsigned Dimension>
template <bool Mutable>
typename OrientedAxes<Dimension>::template Iterator<Mutable>
OrientedAxes<Dimension>::Iterator<Mutable>::operator++(int) noexcept
{
    Iterator copy = *this;
    ++(*this);
    return copy;
}


template <unsigned Dimension>
template <bool Mutable>
typename OrientedAxes<Dimension>::template Iterator<Mutable>&
OrientedAxes<Dimension>::Iterator<Mutable>::operator+=(difference_type offset) noexcept
{
    _index += offset;
    return *this;
}


template <unsigned Dimension>
template <bool Mutable>
typename OrientedAxes<Dimension>::template Iterator<Mutable>
OrientedAxes<Dimension>::Iterator<Mutable>::operator+(difference_type offset) noexcept
{
    Iterator copy = *this;
    copy += offset;
    return copy;
}


template <unsigned Dimension>
template <bool Mutable>
typename OrientedAxes<Dimension>::template Iterator<Mutable>&
OrientedAxes<Dimension>::Iterator<Mutable>::operator--() noexcept
{
    --_index;
    return *this;
}


template <unsigned Dimension>
template <bool Mutable>
typename OrientedAxes<Dimension>::template Iterator<Mutable>
OrientedAxes<Dimension>::Iterator<Mutable>::operator--(int) noexcept
{
    Iterator copy = *this;
    --(*this);
    return copy;
}

template <unsigned Dimension>
template <bool Mutable>
typename OrientedAxes<Dimension>::template Iterator<Mutable>&
OrientedAxes<Dimension>::Iterator<Mutable>::operator-=(difference_type offset) noexcept
{
    _index -= offset;
    return *this;
}


template <unsigned Dimension>
template <bool Mutable>
typename OrientedAxes<Dimension>::template Iterator<Mutable>
OrientedAxes<Dimension>::Iterator<Mutable>::operator-(difference_type offset) noexcept
{
    Iterator copy = *this;
    copy -= offset;
    return copy;
}


template <unsigned Dimension>
template <bool Mutable>
typename OrientedAxes<Dimension>::template Iterator<Mutable>::difference_type
OrientedAxes<Dimension>::Iterator<Mutable>::operator-(Iterator rhs) noexcept
{
    return _index - rhs._index;
}


template <unsigned Dimension>
template <bool Mutable>
template <bool M>
bool OrientedAxes<Dimension>::Iterator<Mutable>::operator==(Iterator<M> right) noexcept
{
    return _pData == right._pData && _index == right._index;
}


template <unsigned Dimension>
template <bool Mutable>
template <bool M>
bool OrientedAxes<Dimension>::Iterator<Mutable>::operator!=(Iterator<M> right) noexcept
{
    return _index != right._index || _pData != right._pData;
}


template <unsigned Dimension>
OrientedAxes<Dimension>::OrientedAxes() noexcept
    : _data(0u)
{
    // By default, the system is not rotated, so all axes align with
    // the original ones.
    for (unsigned iDimension=0u; iDimension<Dimension; ++iDimension) {
        this->at(iDimension) = BoundaryID(iDimension, true);
    }
}


template <unsigned Dimension>
OrientedAxes<Dimension>::OrientedAxes(Ptr<const BoundaryID> itBegin,
                                      size_type size)
{
    CIE_OUT_OF_RANGE_CHECK(size == this->size())
    std::copy(itBegin, itBegin + size, this->begin());
}


template <unsigned Dimension>
OrientedAxes<Dimension>::OrientedAxes(const char axes[2 * Dimension + 1])
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
}


template <unsigned Dimension>
bool OrientedAxes<Dimension>::operator==(OrientedAxes rhs) const noexcept
{
    return this->_data == rhs._data;
}


template <unsigned Dimension>
bool OrientedAxes<Dimension>::operator!=(OrientedAxes rhs) const noexcept
{
    return this->_data != rhs._data;
}


template <unsigned Dimension>
bool OrientedAxes<Dimension>::operator<(OrientedAxes rhs) const noexcept
{
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
typename OrientedAxes<Dimension>::const_reference
OrientedAxes<Dimension>::operator[](size_type index) const
{
    return *(this->begin() + index);
}


template <unsigned Dimension>
typename OrientedAxes<Dimension>::reference
OrientedAxes<Dimension>::operator[](size_type index)
{
    return *(this->begin() + index);
}


template <unsigned Dimension>
typename OrientedAxes<Dimension>::const_reference
OrientedAxes<Dimension>::at(size_type index) const
{
    return *(this->begin() + index);
}


template <unsigned Dimension>
typename OrientedAxes<Dimension>::reference
OrientedAxes<Dimension>::at(size_type index)
{
    return *(this->begin() + index);
}


template <unsigned Dimension>
typename OrientedAxes<Dimension>::const_iterator
OrientedAxes<Dimension>::begin() const noexcept
{
    return const_iterator(_data, 0u);
}


template <unsigned Dimension>
typename OrientedAxes<Dimension>::iterator
OrientedAxes<Dimension>::begin() noexcept
{
    return iterator(_data, 0u);
}


template <unsigned Dimension>
typename OrientedAxes<Dimension>::const_iterator
OrientedAxes<Dimension>::end() const noexcept
{
    return const_iterator(_data, Dimension);
}


template <unsigned Dimension>
typename OrientedAxes<Dimension>::iterator
OrientedAxes<Dimension>::end() noexcept
{
    return iterator(_data, Dimension);
}


template <unsigned Dimension>
void io::GraphML::Serializer<OrientedAxes<Dimension>>::header(Ref<XMLElement> rElement) noexcept
{
    auto descriptionElement = rElement.addChild("desc");

    std::stringstream description;
    description << Dimension << "-dimensional cube spanning [-1, 1]^"
                << Dimension << " in any axis-parallel orientation within "
                << Dimension << "-dimensional space.";
    descriptionElement.setValue(description.view());
}


template <unsigned Dimension>
void io::GraphML::Serializer<OrientedAxes<Dimension>>::operator()(Ref<XMLElement> rElement,
                                                                  Ref<const OrientedAxes<Dimension>> rObject) noexcept
{
    std::stringstream stream;
    stream << rObject;
    rElement.setValue(stream.view());
}


} // namespace cie::fem


#endif
