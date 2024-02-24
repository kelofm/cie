#ifndef CIE_UTILS_COMPILE_TIME_PARAMETER_PACK_SELECT_HPP
#define CIE_UTILS_COMPILE_TIME_PARAMETER_PACK_SELECT_HPP

// --- Utility Includes ---
#include "packages/types/inc/types.hpp"

// --- STL Includes ---
#include <type_traits>
#include <tuple>


namespace cie::ct {


/// @addtogroup cieutils
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


/// @addtogroup cieutils
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


} // namespace cie::ct


#endif