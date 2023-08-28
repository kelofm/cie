#ifndef CIE_UTILS_STL_EXTENSION_FUNCTION_TRAITS_HPP
#define CIE_UTILS_STL_EXTENSION_FUNCTION_TRAITS_HPP

// --- Utility Includes ---
#include "packages/types/inc/types.hpp"

// --- STL Includes ---
#include <functional>
#include <tuple>


namespace cie {


template <class T>
struct FunctionTraits : public std::false_type
{
    static constexpr bool Value = false;
}; // struct FunctionTraits

template <class TReturn, class ...TArguments>
struct FunctionTraits<std::function<TReturn(TArguments...)>> : public std::true_type
{
    static constexpr bool Value = true;
    using Return        = TReturn;
    using result_type   = Return;
    using ArgumentTuple = std::tuple<TArguments...>;

    static constexpr Size NumberOfArguments = sizeof...(TArguments);

    template <Size I_Argument>
    struct Argument
    {
        using Type = typename std::tuple_element<I_Argument, ArgumentTuple>::type;
        using type = Type;
    };
}; // struct FunctionTraits


} // namespace cie


#endif