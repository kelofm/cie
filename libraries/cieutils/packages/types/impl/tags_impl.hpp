#ifndef CIE_UTILS_TAGS_IMPL_HPP
#define CIE_UTILS_TAGS_IMPL_HPP

// --- Utility Includes ---
#include "packages/types/inc/tags.hpp"


namespace cie::tags {


namespace detail {


template <class ...TTags>
constexpr Flags allExcept() noexcept;


template <class TCurrent, class ...TRest>
struct IterateTags
{
    static constexpr Flags allExcept(Ref<Flags> r_flags) noexcept
    {
        r_flags.set(TCurrent::id(), false);
        return IterateTags<TRest...>::allExcept(r_flags);
    }

    template <class TTag>
    static constexpr bool isCompatibleWith() noexcept
    {
        return TCurrent::getCompatibility().test(TTag::id()) && IterateTags<TRest...>::template isCompatibleWith<TTag>();
    }
}; // struct IterateTags


template <class TCurrent>
struct IterateTags<TCurrent>
{
    static constexpr Flags allExcept(Ref<Flags> r_flags) noexcept
    {
        r_flags.set(TCurrent::id(), false);
        return r_flags;
    }

    template <class TTag>
    static constexpr bool isCompatibleWith() noexcept
    {
        return TCurrent::getCompatibility().test(TTag::id());
    }
}; // struct IterateTags


template <class ...TTags>
constexpr Flags allExcept() noexcept
{
    Flags flags;
    flags.set();
    return IterateTags<TTags...>::allExcept(flags);
}


template <std::uint8_t ID, class TSelf>
inline constexpr std::uint8_t
Tag<ID,TSelf>::id() noexcept
{
    return ID;
}


template <std::uint8_t ID, class TSelf>
inline constexpr Flags
Tag<ID,TSelf>::flags() noexcept
{
    Flags flags = 0;
    flags.set(id());
    return flags;
}


template <std::uint8_t ID, class TSelf>
inline constexpr detail::Flags
Tag<ID,TSelf>::getCompatibility() noexcept
{
    Flags flags;
    flags.set();
    return flags;
}


template <std::uint8_t ID, class TSelf>
template <class ...TTags>
inline constexpr bool
Tag<ID,TSelf>::isCompatibleWith() noexcept
{
    return IterateTags<TTags...>::template isCompatibleWith<TSelf>();
}



} // namespace detail


inline constexpr detail::Flags Lazy::getCompatibility() noexcept
{
    return detail::allExcept<Eager>();
}


inline constexpr detail::Flags Eager::getCompatibility() noexcept
{
    return detail::allExcept<Lazy>();
}


inline constexpr detail::Flags Serial::getCompatibility() noexcept
{
    return detail::allExcept<SMP,MPI>();
}


inline constexpr detail::Flags SMP::getCompatibility() noexcept
{
    return detail::allExcept<Serial,MPI>();
}


inline constexpr detail::Flags MPI::getCompatibility() noexcept
{
    return detail::allExcept<Serial,SMP>();
}


inline constexpr detail::Flags Binary::getCompatibility() noexcept
{
    return detail::allExcept<Text>();
}


inline constexpr detail::Flags Text::getCompatibility() noexcept
{
    return detail::allExcept<Binary>();
}


} // namespace cie::tags

#endif