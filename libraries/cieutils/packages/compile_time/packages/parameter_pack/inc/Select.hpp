#ifndef CIE_UTILS_COMPILE_TIME_PARAMETER_PACK_SELECT_HPP
#define CIE_UTILS_COMPILE_TIME_PARAMETER_PACK_SELECT_HPP

// --- Utility Includes ---
#include "packages/types/inc/types.hpp"

// --- STL Includes ---
#include <tuple>
#include <concepts> // std::integral


namespace cie::ct {


/// @ingroup cieutils
template <class ...Ts>
class Select
{
private:
    template <class T>
    struct Tag {using Type = T;};

    using Tuple = std::tuple<Ts...>;

public:
    using Front = typename std::tuple_element<0,Tuple>::type;

    template <Size I>
    requires (I < sizeof...(Ts))
    using At = typename std::tuple_element<I,Tuple>::type;

    using Back = typename decltype((Tag<Ts>{}, ...))::Type;
}; // class Select


/// @ingroup cieutils
template <>
class Select<>
{
public:
    using Front = void;

    // member At is forbidden for empty packs.
    // template <Size I>
    // using At = void;

    using Back = void;
}; // class Select<>


namespace impl {
template <std::integral T, T ...TValues>
struct StripIntegralArray
{}; // struct StripIntegralArray<T,T...>

template <std::integral T, T TLast>
struct StripIntegralArray<T,TLast>
{
    constexpr static T First = TLast;
    using Bottom = std::integer_sequence<T>;
    constexpr static T Last = TLast;
}; // struct StripIntegralArray<T,T>

template <std::integral T, T TFirst, T ...TRest>
struct StripIntegralArray<T,TFirst,TRest...>
{
    constexpr static T First = TFirst;
    using Bottom = std::integer_sequence<T,TRest...>;
    constexpr static T Last = StripIntegralArray<T, TRest...>::Value;
}; // struct StripIntegralArray<T,T,T...>
} // namespace impl


template <std::integral T, T ...Is>
struct SelectIntegral
{
    constexpr static T First = impl::StripIntegralArray<T,Is...>::First;
    using Bottom = typename impl::StripIntegralArray<T,Is...>::Bottom;
    constexpr static T Last = impl::StripIntegralArray<T,Is...>::Last;
}; // struct SelectIntegral


} // namespace cie::ct


#endif