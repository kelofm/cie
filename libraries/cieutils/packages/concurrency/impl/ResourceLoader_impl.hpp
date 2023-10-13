#ifndef CIE_MP_RESOURCE_LOADER_IMPL_HPP
#define CIE_MP_RESOURCE_LOADER_IMPL_HPP

// --- Utility Includes ---
#include "packages/concurrency/inc/ResourceLoader.hpp"
#include "packages/macros/inc/exceptions.hpp"
#include "packages/macros/inc/logging.hpp"


namespace cie::mp {


template <class TResource>
ResourceLoader<TResource>::ResourceLoader(RightRef<Pool> r_pool) noexcept
    : _pool(std::move(r_pool))
{
}


template <class TResource>
ResourceLoader<TResource>::ResourceLoader(Ref<const Pool> r_pool)
    : _pool(r_pool)
{
}


template <class TResource>
template <class ...TArgs>
typename ResourceLoader<TResource>::Value
ResourceLoader<TResource>::load(TArgs&&... r_args)
{
    CIE_BEGIN_EXCEPTION_TRACING

    auto p_promise = std::make_shared<std::promise<TResource>>();
    auto future = p_promise->get_future();

    if (1 < _pool.size()) {
        _pool.queueJob([p_promise = std::move(p_promise), ...r_args = std::forward<TArgs>(r_args)]() mutable {
            CIE_BEGIN_EXCEPTION_TRACING
            CIE_PROFILE_SCOPE
            p_promise->set_value(TResource(std::forward<TArgs>(r_args)...));
            CIE_END_EXCEPTION_TRACING
        });
    } else {
        p_promise->set_value(TResource(std::forward<TArgs>(r_args)...));
    }

    return future;

    CIE_END_EXCEPTION_TRACING
}


} // namespace cie::mp


#endif
