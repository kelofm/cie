#pragma once

// --- STL Includes ---
#include <cstddef>
#include <span>
#include <memory>
#include <atomic>


namespace cie {


struct LinearArena {
    std::span<std::byte> buffer = {};

    std::atomic<std::size_t> state = 0ul;
}; // struct LinearArena


template <class T>
class LinearAllocator {
public:
    using value_type = T;

    using propagate_on_container_move_assignment = std::true_type;

    constexpr LinearAllocator(LinearArena& rArena) noexcept
        :   _pArena(&rArena)
    {}

    template <class U>
    constexpr LinearAllocator(const LinearAllocator<U>& rOther) noexcept
        :   _pArena(rOther._pArena)
    {}

    [[nodiscard]] T* allocate(std::size_t count) noexcept {
        if (_pArena->buffer.size() < _pArena->state) std::terminate();
        const std::size_t byteCount = count * sizeof(T);
        void* itState = _pArena->buffer.data() + _pArena->state;
        std::size_t availableByteCount = static_cast<std::size_t>(_pArena->buffer.size() - _pArena->state);

        void* itNewState = std::align(
            alignof(T),
            byteCount,
            itState,
            availableByteCount);
        if (itNewState == nullptr)
            std::terminate();

        _pArena->state = static_cast<std::size_t>(static_cast<std::byte*>(itNewState) - _pArena->buffer.data());
        return static_cast<T*>(itNewState);
    }

    constexpr void deallocate(
        T*,
        std::size_t) noexcept
    {}

    template <class U>
    [[nodiscard]] constexpr bool operator==(const LinearAllocator<U>& rOther) const noexcept {
        return _pArena->buffer.data() == rOther._pArena->buffer.data() && _pArena->buffer.size() == rOther.arena().size();
    }

    template <class U>
    [[nodiscard]] constexpr bool operator!=(const LinearAllocator<U>& rOther) const noexcept {
        return !((*this) == rOther);
    }

private:
    template <class U>
    friend class LinearAllocator;

    LinearAllocator() = delete;

    LinearArena* _pArena;
}; // class LinearAllocator


} // namespace cie
