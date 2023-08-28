#ifndef CIE_CIEUTILS_TYPES_IMPL_HPP
#define CIE_CIEUTILS_TYPES_IMPL_HPP

// --- Utility Includes ---
#include "packages/types/inc/types.hpp"


namespace cie {


template <class T>
Ptr<const T> makePtrTo(Ref<const T> r_t)
{
    return &r_t;
}


template <class T>
Ptr<T> makePtrTo(Ref<T> r_t)
{
    return &r_t;
}


} // namespace cie


#endif