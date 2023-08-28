// --- Internal Includes ---
#include "packages/observer/inc/Subject.hpp"
#include "packages/macros/inc/exceptions.hpp"

// --- STL Includes ---
#include <memory>
#include <vector>


namespace cie::utils::observer {


Subject::Subject(const SubjectPtr&)
{
}


SubjectPtr Subject::operator=(const SubjectPtr&)
{
    return shared_from_this();
}


Subject::~Subject()
{
    for (auto it=_observers.begin(); it!=_observers.end(); ++it)
        it->reset();
}


void Subject::notifyObservers()
{
    CIE_BEGIN_EXCEPTION_TRACING

    for (auto it=_observers.begin(); it!=_observers.end(); ++it)
        (*it)->onSubjectChange();

    CIE_END_EXCEPTION_TRACING
}


int Subject::attachObserver( ObserverPtr observer )
{
    CIE_BEGIN_EXCEPTION_TRACING

    int id = observerID(observer);

    if ( id == -1 )
    {
        _observers.push_back(observer);
        id = _observers.size()-1;
    }

    return id;

    CIE_END_EXCEPTION_TRACING
}


void Subject::detachObserver( ObserverPtr observer )
{
    CIE_BEGIN_EXCEPTION_TRACING

    int id = observerID(observer);

    if ( observer->subject().get() == this )
        observer->subject().reset();

    if ( id != -1 )
        _observers.erase( _observers.begin()+id );

    CIE_END_EXCEPTION_TRACING
}


int Subject::observerID( ObserverPtr observer )
{
    CIE_BEGIN_EXCEPTION_TRACING

    int id = -1;

    for ( size_t i=0; i<_observers.size(); ++i )
        if ( _observers[i].get() == observer.get() )
        {
            id = (int)i;
            break;
        }

    return id;

    CIE_END_EXCEPTION_TRACING
}


const std::vector<ObserverPtr>& Subject::observers() const
{
    return _observers;
}


} // namespace cie::utils::observer