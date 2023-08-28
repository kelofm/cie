#ifndef CIE_UTILS_PARAMETER_PACK_UNIQUE_HPP
#define CIE_UTILS_PARAMETER_PACK_UNIQUE_HPP

// --- Utility Includes ---
#include "packages/compile_time/packages/parameter_pack/inc/Range.hpp"
#include "packages/compile_time/packages/parameter_pack/inc/Match.hpp"
#include "packages/compile_time/packages/parameter_pack/inc/Size.hpp"


namespace cie::ct {


///@addtogroup cieutils
///@{

namespace detail {
template <class T>
struct IsUnique {};

template <class ...TArgs>
struct IsUnique<std::tuple<TArgs...>>
{
    inline static constexpr bool isUnique()
    requires (!Empty<TArgs...>)
    {
        using Current = typename Select<TArgs...>::Front;
        using Rest = typename Range<TArgs...>::template Subrange<1,PackSize<TArgs...>>::Tuple;
        constexpr const bool current = MatchTuple<std::tuple<Current>>::template None<Rest>;
        return current ? IsUnique<Rest>::isUnique() : false;
    }

    inline static constexpr bool isUnique()
    requires Empty<TArgs...>
    {return true;}
}; // class IsUnique<tuple>
} // namespace detail


/// True if the provided parameter pack has no duplicate types.
/// @note The parameter pack must not contain @a void.
/// @note O^2 implementation :'(
template <class ...TArgs>
inline constexpr const bool IsUnique = detail::IsUnique<std::tuple<TArgs...>>::isUnique();

///@}


} // namespace cie::ct


#endif