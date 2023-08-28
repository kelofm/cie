#ifndef CIE_FEM_MATHS_CARTESIAN_PRODUCT_IMPL_HPP
#define CIE_FEM_MATHS_CARTESIAN_PRODUCT_IMPL_HPP

// --- FEM Includes ---
#include "packages/maths/inc/CartesianProduct.hpp"

// --- Utility Includes ---
#include "packages/macros/inc/checks.hpp"

// --- STL Includes ---
#include <algorithm>


namespace cie::fem::maths {


template <unsigned Dimension>
template <concepts::Integer TValue>
inline bool
CartesianProduct<Dimension>::next(unsigned numberOfStates,
                                  Ptr<TValue> it_stateBegin)
{
    CIE_OUT_OF_RANGE_CHECK(numberOfStates != 0)
    constexpr TValue lowerBound       = 0;
    const TValue upperBound           = numberOfStates - 1;
    Ptr<const TValue> it_stateEnd = it_stateBegin + Dimension;

    for (Ptr<TValue> it=it_stateBegin; it!=it_stateEnd; ++it) {
        if (*it != upperBound) { // <== found a component to increment
            ++(*it);
            return true;
        } else { // overflow: reset every component until this one
            *it = lowerBound;
        }
    }

    // Could not find component to increment
    return false;
}


} // namespace cie::fem::maths


#endif
