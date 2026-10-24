#pragma once

// --- STL Includes ---
#include <cstddef>
#include <cstdint>
#include <span>


namespace cie {


struct LinearArena {
    std::span<std::byte> buffer = {};

    /// @brief Offset of the next allocation.
    /// @details Deliberately not atomic: an arena is not meant to be shared
    ///          between threads or work items. Each consumer (e.g. a SYCL work
    ///          item deserializing an object) owns its arena exclusively.
    std::size_t state = 0ul;
}; // struct LinearArena


template <class T>
class LinearAllocator {
public:
    using value_type = T;

    using propagate_on_container_move_assignment = std::true_type;

    /// @brief Construct an allocator holding no arena.
    /// @details Allocating through a default-constructed allocator terminates
    ///          the program. Such allocators exist in default-constructed
    ///          containers that are reassigned before use (e.g. during
    ///          deserialization), which makes them usable in SYCL kernels.
    constexpr LinearAllocator() noexcept = default;

    constexpr LinearAllocator(LinearArena& rArena) noexcept
        :   _pArena(&rArena)
    {}

    template <class U>
    constexpr LinearAllocator(const LinearAllocator<U>& rOther) noexcept
        :   _pArena(rOther._pArena)
    {}

    [[nodiscard]] T* allocate(std::size_t count) noexcept {
        if (_pArena == nullptr || _pArena->buffer.size() < _pArena->state)
            __builtin_trap();
        const std::size_t byteCount = count * sizeof(T);
        const std::uintptr_t currentAddress = reinterpret_cast<std::uintptr_t>(
            _pArena->buffer.data() + _pArena->state);
        const std::uintptr_t alignedAddress = (currentAddress + alignof(T) - 1ul)
            & ~(static_cast<std::uintptr_t>(alignof(T)) - 1ul);
        const std::size_t alignmentPadding = static_cast<std::size_t>(alignedAddress - currentAddress);
        if (alignmentPadding > _pArena->buffer.size() - _pArena->state
            || byteCount > _pArena->buffer.size() - _pArena->state - alignmentPadding)
            __builtin_trap();

        _pArena->state += alignmentPadding + byteCount;
        return reinterpret_cast<T*>(alignedAddress);
    }

    constexpr void deallocate(
        T*,
        std::size_t) noexcept
    {}

    template <class U>
    [[nodiscard]] constexpr bool operator==(const LinearAllocator<U>& rOther) const noexcept {
        return _pArena == rOther._pArena;
    }

    template <class U>
    [[nodiscard]] constexpr bool operator!=(const LinearAllocator<U>& rOther) const noexcept {
        return !((*this) == rOther);
    }

private:
    template <class U>
    friend class LinearAllocator;

    LinearArena* _pArena = nullptr;
}; // class LinearAllocator


} // namespace cie
