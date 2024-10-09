#ifndef CIE_MP_RESOURCE_LOADER_HPP
#define CIE_MP_RESOURCE_LOADER_HPP

// --- Internal Includes ---
#include "packages/concurrency/inc/ThreadPool.hpp"

// --- STL Includes ---
#include <future>


namespace cie::mp {


/// @ingroup cieutils
template <class TResource>
class ResourceLoader
{
public:
    using Pool = ThreadPool<>;

    using Value = std::future<TResource>;

public:
    ResourceLoader(RightRef<Pool> r_pool) noexcept;

    ResourceLoader(Ref<const Pool> r_pool);

    ResourceLoader(ResourceLoader&& r_rhs) noexcept = default;

    ResourceLoader(const ResourceLoader& r_rhs) = delete;

    ResourceLoader& operator=(ResourceLoader&& r_rhs) = default;

    ResourceLoader& operator=(const ResourceLoader& r_rhs) = delete;

    template <class ...TArgs>
    [[nodiscard]] Value load(TArgs&&... r_args);

private:
    Pool _pool;
}; // class ResourceLoader


} // namespace cie::mp

#include "packages/concurrency/impl/ResourceLoader_impl.hpp"

#endif