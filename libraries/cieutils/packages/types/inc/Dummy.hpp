#pragma once


namespace cie {


struct Dummy {
    constexpr Dummy() noexcept = default;

    template <class T>
    constexpr Dummy(T&&) noexcept {}
}; // struct Dummy


} // namespace cie
