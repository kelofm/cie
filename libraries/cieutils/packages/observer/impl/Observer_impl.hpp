#ifndef CIEUTILS_ABSOBSERVER_IMPL_HPP
#define CIEUTILS_ABSOBSERVER_IMPL_HPP

// --- Internal Includes ---
#include "packages/macros/inc/exceptions.hpp"


namespace cie::utils::observer {


template <class DerivedObserver>
std::shared_ptr<DerivedObserver> Observer::as()
{
    CIE_BEGIN_EXCEPTION_TRACING

    return std::dynamic_pointer_cast<DerivedObserver>( shared_from_this() );

    CIE_END_EXCEPTION_TRACING
}


} // namespace cie::utils::observer

#endif