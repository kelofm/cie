// --- Internal Ipmorts ---
#include "packages/observer/inc/Observer.hpp"
#include "packages/observer/inc/Subject.hpp"
#include "packages/macros/inc/exceptions.hpp"
#include "packages/macros/inc/checks.hpp"

namespace cie::utils::observer {


Observer::Observer() :
    _subject( )
{
}


Observer::Observer( SubjectPtr subject ) :
    _subject( subject )
{
}


Observer::Observer( const ObserverPtr& copy ) :
    _subject( copy->_subject )
{
}


ObserverPtr Observer::operator=( const ObserverPtr& copy )
{
    CIE_BEGIN_EXCEPTION_TRACING

    CIE_CHECK_POINTER( copy )
    _subject = copy->_subject;
    return shared_from_this();

    CIE_END_EXCEPTION_TRACING
}


Observer::~Observer()
{
    _subject.reset();
}


void Observer::detach()
{
    CIE_BEGIN_EXCEPTION_TRACING

    if ( auto p_subject = _subject.lock() )
        p_subject->detachObserver( this->shared_from_this() );

    CIE_END_EXCEPTION_TRACING
}


void Observer::setSubject( SubjectPtr subject )
{
    _subject = subject;
}


SubjectPtr Observer::subject()
{
    return _subject.lock();
}


const SubjectPtr Observer::subject() const
{
    return _subject.lock();
}


} // namespace cie::utils::observer