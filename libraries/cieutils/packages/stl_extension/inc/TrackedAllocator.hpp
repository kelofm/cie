#pragma once

// --- STL Includes ---
#include <memory>
#include <mutex>
#include <atomic>


namespace cie {


class AllocatorStats final {
public:
    [[nodiscard]] std::size_t current() const noexcept {
        return _current;
    }

    [[nodiscard]] std::size_t peak() const noexcept {
        return _peak;
    }

private:
    template <template <class ...> class TAllocator, class T, class ...Ts>
    friend class TrackedAllocator;

    void reportAllocation(std::size_t size) noexcept {
        _current += size;
        std::scoped_lock<std::mutex> lock(_mutex);
        _peak = std::max<std::size_t>(
            _peak,
            _current);
    }

    void reportDeallocation(std::size_t size) noexcept {
        _current -= size;
    }

    std::mutex _mutex;

    std::atomic<std::size_t> _current = 0ul;

    std::atomic<std::size_t> _peak = 0ul;
}; // class AllocatorStats


template <template <class ...> class TAllocator, class T, class ...Ts>
class TrackedAllocator {
public:
    using value_type = T;

    constexpr TrackedAllocator() noexcept
    requires std::is_default_constructible_v<TAllocator<T,Ts...>> = default;

    TrackedAllocator(
        TAllocator<T,Ts...>&& rAllocator,
        AllocatorStats& rStats)
            :   _allocator(std::move(rAllocator)),
                _pStats(&rStats)
    {}

    template <class U, class ...Us>
    TrackedAllocator(const TrackedAllocator<TAllocator,U,Us...> rOther) noexcept
    requires std::is_constructible_v<TAllocator<T,Ts...>,TAllocator<U,Us...>>
        :   _allocator(rOther),
            _pStats(rOther._pStats)
    {}

    T* allocate(std::size_t count) {
        const std::size_t byteCount = count * sizeof(T);
        _pStats->reportAllocation(byteCount);
        return _allocator.allocate(count);
    }

    void deallocate(
        T* pBegin,
        std::size_t count) {
            const std::size_t byteCount = count * sizeof(T);
            _pStats->reportDeallocation(byteCount);
            _allocator.deallocate(pBegin, count);
    }

    template <class U, class ...Us>
    [[nodiscard]] constexpr bool operator==(const TrackedAllocator<TAllocator,U,Us...> rOther) const noexcept {
        return _allocator == rOther._allocator && _pStats == rOther._pStats;
    }

    template <class U, class ...Us>
    [[nodiscard]] constexpr bool operator!=(const TrackedAllocator<TAllocator,U,Us...> rOther) const noexcept {
        return !((*this) == rOther);
    }

private:
    template <template <class ...> class TA, class U, class ...Us>
    friend class TrackedAllocator;

    TAllocator<T,Ts...> _allocator;

    AllocatorStats* _pStats;
}; // class TrackedAllocator


} // namespace cie
