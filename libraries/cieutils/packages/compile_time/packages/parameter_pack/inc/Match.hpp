#ifndef CIE_UTILS_COMPILE_TIME_PARAMETER_PACK_MATCH_HPP
#define CIE_UTILS_COMPILE_TIME_PARAMETER_PACK_MATCH_HPP

// --- Utility Includes ---
#include "packages/types/inc/types.hpp"
#include "packages/compile_time/packages/parameter_pack/inc/Select.hpp"
#include "packages/compile_time/packages/parameter_pack/inc/Size.hpp"

// --- STL Includes ---
#include <type_traits>
#include <tuple>


namespace cie::ct {


///@addtogroup cieutils
///@{


namespace detail {

/// Compare a type to a parameter pack.
template <class T>
class MatchOne
{
private:
    template <Size I, class ...TRefs>
    requires (!Empty<TRefs...>)
    static constexpr bool any()
    {
        if constexpr (std::is_same_v<T,typename Select<TRefs...>::template At<I>>)
            return true;

        if constexpr (I == 0)
            return false;
        else
            return any<I-1,TRefs...>();
    }

    template <Size I, class ...TRefs>
    requires Empty<TRefs...>
    static constexpr bool any() {return false;}

public:
    /// T is in @ref TRefs
    template <class ...TRefs>
    inline static constexpr bool Any = any<PackSize<TRefs...>-1,TRefs...>();

    /// T is not in @ref TRefs
    template <class ...TRefs>
    inline static constexpr bool None = !Any<TRefs...>;
}; // class MatchOne

} // namespace detail


/// Compare two parameter packs.
/// @addtogroup cieutils
template <class ...Ts>
class Match
{
private:
    template <Size I, class ...TRefs>
    requires (!Empty<TRefs...> && !Empty<Ts...>)
    static constexpr bool any()
    {
        using Current = typename Select<Ts...>::template At<I>;
        if constexpr (detail::MatchOne<Current>::template Any<TRefs...>)
            return true;

        if constexpr (I == 0)
            return false;
        else
            return any<I-1,TRefs...>();
    }

    template <Size I, class ...TRefs>
    requires (Empty<TRefs...> || Empty<Ts...>)
    static constexpr bool any() {return false;}

    template <Size I, class ...TRefs>
    requires (!Empty<TRefs...> && !Empty<Ts...>)
    static constexpr bool all()
    {
        using Current = typename Select<Ts...>::template At<I>;
        constexpr const bool current = detail::MatchOne<Current>::template Any<TRefs...>;

        if constexpr (I == 0)
            return current;
        else
            return current ? all<I-1,TRefs...>() : false;
    }

    template <Size I, class ...TRefs>
    requires (Empty<TRefs...> || Empty<Ts...>)
    static constexpr bool all() {return Empty<TRefs...> && Empty<Ts...>;}

    template <Size I, class ...TRefs>
    requires (!Empty<TRefs...> && (PackSize<TRefs...> == PackSize<Ts...>))
    static constexpr bool identical()
    {
        using Current = typename Select<Ts...>::template At<I>;
        constexpr const bool current = std::is_same_v<Current, typename Select<TRefs...>::template At<I>>;

        if constexpr (I == 0)
            return current;
        else
            return current ? identical<I-1,TRefs...>() : false;
    }

    template <Size I, class ...TRefs>
    requires (Empty<TRefs...> || (PackSize<TRefs...> != PackSize<Ts...>))
    static constexpr bool identical() {return PackSize<TRefs...> == PackSize<Ts...>;}

public:
    /// Any of @p Ts is in @p TRefs
    template <class ...TRefs>
    inline static constexpr bool Any = any<PackSize<Ts...>-1,TRefs...>();

    /// None of @p Ts is in @p TRefs
    template <class ...TRefs>
    inline static constexpr bool None = !Any<TRefs...>;

    /// All of @p Ts is present in @p TRefs
    template <class ...TRefs>
    inline static constexpr bool All = all<PackSize<Ts...>-1,TRefs...>();

    /// True if @p Ts and @p TRefs are identical parameter packs (same types in the same order).
    template <class ...TRefs>
    inline static constexpr bool Identical = identical<PackSize<Ts...>-1,TRefs...>();
}; // class Match


/// Compare the parameter packs of two tuples
template <class T>
class MatchTuple {};


/// Compare the parameter packs of two tuples.
/// @addtogroup cieutils
template <class ...Ts>
class MatchTuple<std::tuple<Ts...>>
{
private:
    template <class T>
    struct AnyImpl {};

    template <class ...TTs>
    struct AnyImpl<std::tuple<TTs...>>
    {inline static constexpr bool Value = Match<Ts...>::template Any<TTs...>;};

    template <class T>
    struct AllImpl {};

    template <class ...TTs>
    struct AllImpl<std::tuple<TTs...>>
    {inline static constexpr bool Value = Match<Ts...>::template All<TTs...>;};

    template <class T>
    struct IdenticalImpl {};

    template <class ...TTs>
    struct IdenticalImpl<std::tuple<TTs...>>
    {inline static constexpr bool Value = Match<Ts...>::template Identical<TTs...>;};

public:
    /// True if any type in the tuple is also in the input tuple.
    template <class TTuple>
    inline static constexpr bool Any = AnyImpl<TTuple>::Value;

    /// True if none of the tuple's types are in the input tuple.
    template <class TTuple>
    inline static constexpr bool None = !Any<TTuple>;

    /// True if all of the tuple's types are in the input tuple.
    template <class TTuple>
    inline static constexpr bool All = AllImpl<TTuple>::Value;

    /// True if the two tuples are identical.
    template <class TTuple>
    inline static constexpr bool Identical = IdenticalImpl<TTuple>::Value;
}; // class MatchTuple



template <class TSet, class TExcept>
class Filter {};



template <class ...TSet, class ...TExcept>
class Filter<std::tuple<TSet...>,std::tuple<TExcept...>>
{
private:
    using SetTypes = std::tuple<TSet...>;

    using ExceptTypes = std::tuple<TExcept...>;

public:
    using Type = decltype(std::tuple_cat(
        std::declval<std::conditional_t<
            Match<TSet>::template Any<TExcept...>,
            std::tuple<>,
            std::tuple<TSet>
        >>()...
    ));
}; // class Filter


///@}


} // namespace cie::ct


#endif