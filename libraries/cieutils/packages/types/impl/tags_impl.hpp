#ifndef CIE_UTILS_TAGS_IMPL_HPP
#define CIE_UTILS_TAGS_IMPL_HPP

// --- Utility Includes ---
#include "packages/types/inc/tags.hpp"


namespace cie::tags {


namespace detail {


template <class ...TTags>
constexpr Flags allExcept();


template <class TCurrent, class ...TRest>
struct IterateTags
{
    static constexpr Flags allExcept(Ref<Flags> r_flags)
    {
        r_flags.set(TCurrent::id(), false);
        return IterateTags<TRest...>::allExcept(r_flags);
    }

    template <class TTag>
    static constexpr bool isCompatibleWith()
    {
        return TCurrent::getCompatibility().test(TTag::id()) && IterateTags<TRest...>::template isCompatibleWith<TTag>();
    }
}; // struct IterateTags


template <class TCurrent>
struct IterateTags<TCurrent>
{
    static constexpr Flags allExcept(Ref<Flags> r_flags)
    {
        r_flags.set(TCurrent::id(), false);
        return r_flags;
    }

    template <class TTag>
    static constexpr bool isCompatibleWith()
    {
        return TCurrent::getCompatibility().test(TTag::id());
    }
}; // struct IterateTags


template <class ...TTags>
constexpr Flags allExcept()
{
    Flags flags;
    flags.set();
    return IterateTags<TTags...>::allExcept(flags);
}


template <Size ID, class TSelf>
inline constexpr Size
Tag<ID,TSelf>::id()
{
    return ID;
}


template <Size ID, class TSelf>
inline constexpr detail::Flags
Tag<ID,TSelf>::getCompatibility()
{
    Flags flags;
    flags.set();
    return flags;
}


template <Size ID, class TSelf>
template <class ...TTags>
inline constexpr bool
Tag<ID,TSelf>::isCompatibleWith()
{
    return IterateTags<TTags...>::template isCompatibleWith<TSelf>();
}



} // namespace detail


inline constexpr detail::Flags Lazy::getCompatibility()
{
    return detail::allExcept<Eager>();
}


inline constexpr detail::Flags Eager::getCompatibility()
{
    return detail::allExcept<Lazy>();
}


inline constexpr detail::Flags Serial::getCompatibility()
{
    return detail::allExcept<SMP,MPI>();
}


inline constexpr detail::Flags SMP::getCompatibility()
{
    return detail::allExcept<Serial,MPI>();
}


inline constexpr detail::Flags MPI::getCompatibility()
{
    return detail::allExcept<Serial,SMP>();
}


inline constexpr detail::Flags Binary::getCompatibility()
{
    return detail::allExcept<Text>();
}


inline constexpr detail::Flags Text::getCompatibility()
{
    return detail::allExcept<Binary>();
}


} // namespace cie::tags

#endif