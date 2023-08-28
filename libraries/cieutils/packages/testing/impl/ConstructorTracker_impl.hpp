#ifndef CIE_UTILS_TESTING_CONSTRUCTOR_TRACKER_IMPL_HPP
#define CIE_UTILS_TESTING_CONSTRUCTOR_TRACKER_IMPL_HPP

// --- Internal Includes ---
#include "packages/macros/inc/exceptions.hpp"


namespace cie::utils {


template <class Tag>
Size ConstructorTracker<Tag>::_numberOfActiveInstances = 0;

template <class Tag>
Size ConstructorTracker<Tag>::_numberOfDefaultConstructorCalls = 0;

template <class Tag>
Size ConstructorTracker<Tag>::_numberOfMoveConstructorCalls = 0;

template <class Tag>
Size ConstructorTracker<Tag>::_numberOfCopyConstructorCalls = 0;

template <class Tag>
Size ConstructorTracker<Tag>::_numberOfMoveAssignmentCalls = 0;

template <class Tag>
Size ConstructorTracker<Tag>::_numberOfCopyAssignmentCalls = 0;

template <class Tag>
Size ConstructorTracker<Tag>::_numberOfDestructorCalls = 0;

template <class Tag>
std::mutex ConstructorTracker<Tag>::_mutex = std::mutex();


template <class Tag>
ConstructorTracker<Tag>::ConstructorTracker()
{
    std::scoped_lock<std::mutex> lock(_mutex);
    ConstructorTracker<Tag>::_numberOfActiveInstances++;
    ConstructorTracker<Tag>::_numberOfDefaultConstructorCalls++;
}


template <class Tag>
ConstructorTracker<Tag>::ConstructorTracker(ConstructorTracker&&)
{
    std::scoped_lock<std::mutex> lock(_mutex);
    ConstructorTracker<Tag>::_numberOfMoveConstructorCalls++;
    ConstructorTracker<Tag>::_numberOfActiveInstances++;
}


template <class Tag>
ConstructorTracker<Tag>::ConstructorTracker(const ConstructorTracker&)
{
    std::scoped_lock<std::mutex> lock(_mutex);
    ConstructorTracker<Tag>::_numberOfCopyConstructorCalls++;
    ConstructorTracker<Tag>::_numberOfActiveInstances++;
}


template <class Tag>
ConstructorTracker<Tag>& ConstructorTracker<Tag>::operator=(ConstructorTracker&&)
{
    std::scoped_lock<std::mutex> lock(_mutex);
    ConstructorTracker<Tag>::_numberOfMoveAssignmentCalls++;
    return *this;
}


template <class Tag>
ConstructorTracker<Tag>& ConstructorTracker<Tag>::operator=(const ConstructorTracker&)
{
    std::scoped_lock<std::mutex> lock(_mutex);
    ConstructorTracker<Tag>::_numberOfCopyAssignmentCalls++;
    return *this;
}


template <class Tag>
ConstructorTracker<Tag>::~ConstructorTracker()
{
    std::scoped_lock<std::mutex> lock(_mutex);
    ConstructorTracker<Tag>::_numberOfActiveInstances--;
    ConstructorTracker<Tag>::_numberOfDestructorCalls++;
}


template <class Tag>
Size ConstructorTracker<Tag>::numberOfActiveInstances()
{
    std::scoped_lock<std::mutex> lock(_mutex);
    return ConstructorTracker<Tag>::_numberOfActiveInstances;
}


template <class Tag>
Size ConstructorTracker<Tag>::numberOfDefaultConstructorCalls()
{
    std::scoped_lock<std::mutex> lock(_mutex);
    return ConstructorTracker<Tag>::_numberOfDefaultConstructorCalls;
}


template <class Tag>
Size ConstructorTracker<Tag>::numberOfMoveConstructorCalls()
{
    std::scoped_lock<std::mutex> lock(_mutex);
    return ConstructorTracker<Tag>::_numberOfMoveConstructorCalls;
}


template <class Tag>
Size ConstructorTracker<Tag>::numberOfCopyConstructorCalls()
{
    std::scoped_lock<std::mutex> lock(_mutex);
    return ConstructorTracker<Tag>::_numberOfCopyConstructorCalls;
}


template <class Tag>
Size ConstructorTracker<Tag>::numberOfMoveAssignmentCalls()
{
    std::scoped_lock<std::mutex> lock(_mutex);
    return ConstructorTracker<Tag>::_numberOfMoveAssignmentCalls;
}


template <class Tag>
Size ConstructorTracker<Tag>::numberOfCopyAssignmentCalls()
{
    std::scoped_lock<std::mutex> lock(_mutex);
    return ConstructorTracker<Tag>::_numberOfCopyAssignmentCalls;
}


template <class Tag>
Size ConstructorTracker<Tag>::numberOfDestructorCalls()
{
    std::scoped_lock<std::mutex> lock(_mutex);
    return ConstructorTracker<Tag>::_numberOfDestructorCalls;
}


template <class Tag>
void ConstructorTracker<Tag>::reset()
{
    std::scoped_lock<std::mutex> lock(_mutex);
    if (ConstructorTracker<Tag>::_numberOfActiveInstances != 0)
        CIE_THROW(Exception, "Cannot reset counters while instances exist")

    ConstructorTracker<Tag>::_numberOfActiveInstances = 0;
    ConstructorTracker<Tag>::_numberOfDefaultConstructorCalls = 0;
    ConstructorTracker<Tag>::_numberOfMoveConstructorCalls = 0;
    ConstructorTracker<Tag>::_numberOfCopyConstructorCalls = 0;
    ConstructorTracker<Tag>::_numberOfMoveAssignmentCalls = 0;
    ConstructorTracker<Tag>::_numberOfCopyAssignmentCalls = 0;
    ConstructorTracker<Tag>::_numberOfDestructorCalls = 0;
}


} // namespace cie::utils

#endif