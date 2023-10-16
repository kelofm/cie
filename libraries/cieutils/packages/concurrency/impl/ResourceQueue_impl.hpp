#ifndef CIE_MP_RESOURCE_QUEUE_IMPL_HPP
#define CIE_MP_RESOURCE_QUEUE_IMPL_HPP

// --- Utility Includes ---
#include "packages/concurrency/inc/ResourceQueue.hpp"
#include "packages/macros/inc/exceptions.hpp"
#include "packages/macros/inc/checks.hpp"


namespace cie::mp {


template <class TResource>
ResourceQueue<TResource>::ResourceQueue(Ref<ThreadPoolBase> r_pool)
    : _loader(r_pool),
      _queue(),
      _mutex()
{
}


template <class TResource>
template <class ...TArgs>
void ResourceQueue<TResource>::emplace(TArgs&&... r_args)
{
    Lock lock(_mutex);
    _queue.push(_loader.load(std::forward<TArgs>(r_args)...));
}


template <class TResource>
TResource
ResourceQueue<TResource>::pop()
{
    Lock lock(_mutex);
    CIE_OUT_OF_RANGE_CHECK(!this->empty(lock))
    TResource value(std::move(_queue.front().get()));
    _queue.pop();
    return value;
}


template <class TResource>
Ref<const TResource>
ResourceQueue<TResource>::back() const
{
    Lock lock(_mutex);
    CIE_OUT_OF_RANGE_CHECK(!this->empty(lock))
    return _queue.front();
}


template <class TResource>
bool ResourceQueue<TResource>::empty() const
{
    return this->empty(Lock(_mutex));
}


template <class TResource>
bool ResourceQueue<TResource>::empty(Ref<const Lock>) const
{
    return _queue.empty();
}


template <class TResource>
Size ResourceQueue<TResource>::size() const
{
    return this->size(Lock(_mutex));
}


template <class TResource>
Size ResourceQueue<TResource>::size(Ref<const Lock>) const
{
    return _queue.size();
}


} // namespace cie::mp


#endif