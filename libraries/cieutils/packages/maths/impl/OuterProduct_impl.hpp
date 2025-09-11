#ifndef CIE_FEM_MATHS_CARTESIAN_PRODUCT_IMPL_HPP
#define CIE_FEM_MATHS_CARTESIAN_PRODUCT_IMPL_HPP

// --- FEM Includes ---
#include "packages/maths/inc/OuterProduct.hpp"

// --- Utility Includes ---
#include "packages/macros/inc/checks.hpp"

// --- STL Includes ---
#include <algorithm>


namespace cie::maths {


template <unsigned Dimension>
template <concepts::Integer TValue>
inline bool
OuterProduct<Dimension>::next(unsigned numberOfStates,
                              Ptr<TValue> itStateBegin)
{
    CIE_OUT_OF_RANGE_CHECK(0 < numberOfStates)
    return OuterProduct<Dimension>::next(static_cast<TValue>(0),
                                         static_cast<TValue>(numberOfStates - 1u),
                                         itStateBegin);
}


template <unsigned Dimension>
template <std::incrementable TValue>
bool OuterProduct<Dimension>::next(Ref<const TValue> firstState,
                                   Ref<const TValue> lastState,
                                   Ptr<TValue> itStateBegin)
{
    Ptr<const TValue> itStateEnd = itStateBegin + Dimension;

    for (Ptr<TValue> it=itStateBegin; it!=itStateEnd; ++it) {
        if (*it != lastState) { // <== found a component to increment
            ++(*it);
            return true;
        } else { // overflow: reset every component until the next incrementable is found
            *it = firstState;
        }
    }

    // Could not find component to increment
    return false;
}


template <unsigned Dimension>
template <std::incrementable TValue,
          concepts::Iterator<TValue> TItBegin,
          concepts::Iterator<TValue> TItEnd>
bool OuterProduct<Dimension>::next(TItBegin itOptionsBegin,
                                   TItEnd itOptionsEnd,
                                   Ptr<TValue> itStateBegin)
{
    if (itOptionsBegin == itOptionsEnd) return false;

    Ref<const TValue> rLastOption = *(itOptionsEnd - 1);
    Ptr<const TValue> itStateEnd = itStateBegin + Dimension;

    for (Ptr<TValue> it=itStateBegin; it!=itStateEnd; ++it) {
        if (*it != rLastOption) { // <== found a component to increment
            const auto itOption = std::find(itOptionsBegin,
                                            itOptionsEnd,
                                            *it);
            *it = itOption[1];
            return true;
        } else { // overflow: reset every component until the next incrementable is found
            *it = *itOptionsBegin;
        }
    }

    // Could not find component to increment
    return false;
}


template <unsigned Dimension>
template <std::incrementable TValue>
bool OuterProduct<Dimension>::next(Ptr<const Ptr<const TValue>> itOptionBeginsArray,
                                   Ptr<const Ptr<const TValue>> itOptionEndsArray,
                                   Ptr<TValue> itStateBegin)
{
    for (unsigned iComponent=0u; iComponent<Dimension; ++iComponent) {
        Ptr<const TValue> itOptionsBegin = itOptionBeginsArray[iComponent];
        Ptr<const TValue> itOptionsEnd = itOptionEndsArray[iComponent];

        if (itOptionsBegin == itOptionsEnd) {
            return false;
        }

        Ref<const TValue> rLastOption = *(itOptionsEnd - 1);
        Ptr<TValue> it = itStateBegin + iComponent;

        if (*it != rLastOption) { // <== found a component to increment
            const auto itOption = std::find(itOptionsBegin,
                                            itOptionsEnd,
                                            *it);
            *it = itOption[1];
            return true;
        } else { // overflow: reset every component until the next incrementable is found
            *it = *itOptionsBegin;
        }
    }

    // Could not find component to increment
    return false;
}


} // namespace cie::maths


#endif
