#ifndef CIE_MP_RESOURCE_QUEUE_HPP
#define CIE_MP_RESOURCE_QUEUE_HPP

// --- Utility Includes ---
#include "packages/concurrency/inc/ResourceLoader.hpp"
#include "packages/concurrency/inc/ThreadPoolBase.hpp"

// --- STL Includes ---
#include <queue>
#include <mutex>


namespace cie::mp {


/// @ingroup cieutils
template <class TResource>
class ResourceQueue
{
public:
    using Loader = ResourceLoader<TResource>;

public:
    ResourceQueue(Ref<ThreadPoolBase> r_pool);

    ResourceQueue(ResourceQueue&& r_rhs) = default;

    ResourceQueue(const ResourceQueue& r_rhs) = delete;

    ResourceQueue& operator=(ResourceQueue&& r_rhs) = default;

    ResourceQueue& operator=(const ResourceQueue& r_rhs) = delete;

    template <class ...TArgs>
    void emplace(TArgs&&... r_args);

    TResource pop();

    Ref<const TResource> back() const;

    bool empty() const;

    Size size() const;

private:
    using Mutex = std::mutex;

    using Lock = std::scoped_lock<Mutex>;

    bool empty(Ref<const Lock>) const;

    Size size(Ref<const Lock>) const;

private:
    Loader _loader;

    std::queue<typename Loader::Value> _queue;

    mutable Mutex _mutex;
}; // class ResourceQueue


} // namespace cie::mp

#include "packages/concurrency/impl/ResourceQueue_impl.hpp"

#endif