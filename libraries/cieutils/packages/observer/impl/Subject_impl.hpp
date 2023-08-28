#ifndef CIEUTILS_ABSSUBJECT_IMPL_HPP
#define CIEUTILS_ABSSUBJECT_IMPL_HPP

// --- Internal Includes ---
#include "packages/macros/inc/exceptions.hpp"


namespace cie::utils::observer {


template <class DerivedSubject>
std::shared_ptr<DerivedSubject> Subject::as()
{
    CIE_BEGIN_EXCEPTION_TRACING

    return std::dynamic_pointer_cast<DerivedSubject>( this->shared_from_this() );

    CIE_END_EXCEPTION_TRACING
}


} // namespace cie::utils::observer

#endif