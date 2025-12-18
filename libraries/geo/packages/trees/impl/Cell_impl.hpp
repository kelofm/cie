#ifndef CIE_GEO_TREES_CELL_IMPL_HPP
#define CIE_GEO_TREES_CELL_IMPL_HPP

// --- Utility Includes ---
#include "packages/exceptions/inc/exception.hpp"
#include "packages/stl_extension/inc/resize.hpp"
#include "packages/maths/inc/power.hpp"
#include "packages/macros/inc/checks.hpp"


namespace cie::geo {


/* --- Cube cell --- */

// 2 children per dimension
template <concepts::Cube TPrimitive>
const typename Cell<TPrimitive>::IndexConverter
    Cell<TPrimitive>::_childIndexConverter(2);


template <concepts::Cube TPrimitive>
Cell<TPrimitive>::Cell(const typename Cell<TPrimitive>::Point& base,
                       typename Cell<TPrimitive>::Coordinate length) :
    Cell<TPrimitive>::CellBase(base, length)
{
}


template <concepts::Cube TPrimitive>
inline typename Cell<TPrimitive>::primitive_constructor_container
Cell<TPrimitive>::split_internal(const typename Cell<TPrimitive>::Point&)
{
    CIE_BEGIN_EXCEPTION_TRACING

    constexpr Size numberOfChildren = intPow(Size(2),Cell<TPrimitive>::Dimension);
    typename Cell<TPrimitive>::primitive_constructor_container constructorArgumentsContainer;
    utils::reserve(constructorArgumentsContainer, numberOfChildren);
    typename Cell<TPrimitive>::Point tempBase;

    for (Size iChild=0; iChild < numberOfChildren; ++iChild) {
        for (Size iDimension=0; iDimension<Cell<TPrimitive>::Dimension; ++iDimension) {
            if (_childIndexConverter.convert(iChild)[iDimension] == 0)
                tempBase[iDimension]   = this->base()[iDimension];
            else
                tempBase[iDimension]   = this->base()[iDimension] + this->_length/2.0;
        }

        constructorArgumentsContainer.emplace_back(tempBase, this->_length/2.0);
    }

    return constructorArgumentsContainer;

    CIE_END_EXCEPTION_TRACING
}



template <concepts::Cube TPrimitive>
inline typename Cell<TPrimitive>::primitive_constructor_container
Cell<TPrimitive>::split(const typename Cell<TPrimitive>::Point& r_point)
{
    return AbsCell<TPrimitive>::split(r_point);
}


template <concepts::Cube TPrimitive>
inline typename Cell<TPrimitive>::primitive_constructor_container
Cell<TPrimitive>::split()
{
    return this->split_internal(typename Cell<TPrimitive>::Point());
}


/* --- Box cell --- */

// 2 children per dimension
template <concepts::Box TPrimitive>
const typename Cell<TPrimitive>::IndexConverter
    Cell<TPrimitive>::_childIndexConverter(2);


template <concepts::Box TPrimitive>
Cell<TPrimitive>::Cell(const typename Cell<TPrimitive>::Point& base,
                       const typename Cell<TPrimitive>::Point& lengths) :
    Cell<TPrimitive>::CellBase(base, lengths)
{
}


template <concepts::Box TPrimitive>
inline typename Cell<TPrimitive>::primitive_constructor_container
Cell<TPrimitive>::split_internal(const typename Cell<TPrimitive>::Point& point)
{
    typename Cell<TPrimitive>::primitive_constructor_container constructorArgumentsContainer;

    constexpr const Size numberOfChildren = intPow(Size(2),Cell<TPrimitive>::Dimension);
    utils::reserve(constructorArgumentsContainer, numberOfChildren);

    typename Cell<TPrimitive>::Point tempBase, tempLengths;

    for (Size iChild=0; iChild < numberOfChildren; ++iChild) {
        for (Size iDimension=0; iDimension<Cell<TPrimitive>::Dimension; ++iDimension) {
            if (_childIndexConverter.convert(iChild)[iDimension] == 0) {
                tempLengths[iDimension]  = point[iDimension] - this->base()[iDimension];
                tempBase[iDimension]     = this->base()[iDimension];
            } else {
                tempLengths[iDimension]  = (this->base()[iDimension] + this->lengths()[iDimension]) - point[iDimension];
                tempBase[iDimension]     = point[iDimension];
            }
        }

        constructorArgumentsContainer.emplace_back(tempBase, tempLengths);
    }

    return constructorArgumentsContainer;
}


} // namespace cie::geo


#endif