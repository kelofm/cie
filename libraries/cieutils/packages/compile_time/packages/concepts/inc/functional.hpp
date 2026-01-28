#ifndef CIE_UTILS_CONCEPTS_FUNCTIONAL_HPP
#define CIE_UTILS_CONCEPTS_FUNCTIONAL_HPP

// --- Utility Includes ---
#include "packages/stl_extension/inc/FunctionTraits.hpp"

// --- STL Includes ---
#include <concepts>
#include <functional>
#include <type_traits>
#include <tuple>


namespace cie::concepts {


template <class T>
concept STLFunction
= FunctionTraits<T>::value;


template <class T>
concept RawFunction
= std::is_function<T>::value;


template <class T>
concept Function
= STLFunction<T> || RawFunction<T>;


template <class F, class A>
concept FunctionWithArgument
=  STLFunction<F>
&& FunctionTraits<F>::NumberOfArguments == 1
&& std::is_same_v<A, typename FunctionTraits<F>::template Argument<0>::type>;

namespace detail {
template <Function TFunction, class ...TArguments>
class IsFunctionWithArguments
{
private:
    template <Size I, bool All>
    struct IsFunctionWithArgumentAtRecursive
    {
        static constexpr bool Value = All && std::is_same_v<
            typename std::tuple_element<I,std::tuple<TArguments...>>::type,
            typename std::tuple_element<I,typename FunctionTraits<TFunction>::ArgumentTuple>::type
        >;

        static constexpr bool value() {return IsFunctionWithArgumentAtRecursive<I-1,Value>::value();}
    }; // struct IsFunctionWithArgumentAtRecursive

    template <bool All>
    struct IsFunctionWithArgumentAtRecursive<0,All>
    {
        static constexpr bool Value = All && std::is_same_v<
            typename std::tuple_element<0,std::tuple<TArguments...>>::type,
            typename std::tuple_element<0,typename FunctionTraits<TFunction>::ArgumentTuple>::type
        >;

        static constexpr bool value() {return Value;}
    }; // struct isFunctionWithArgumentAtRecursive specialized for 0

    using Traits = FunctionTraits<TFunction>;

public:
    static constexpr bool Value = IsFunctionWithArgumentAtRecursive<Traits::NumberOfArguments-1,true>::value();
}; // struct IsFunctionWithArguments
} // namespace detail


template <class TFunction, class ...TArguments>
concept FunctionWithArguments
=  Function<TFunction>
&& FunctionTraits<TFunction>::NumberOfArguments == sizeof...(TArguments)
&& detail::IsFunctionWithArguments<TFunction,TArguments...>::Value;


template <class TFunction, class ...TArgs>
concept CallableWith = requires (TFunction function, TArgs... arguments) {
    {function(std::forward<TArgs>(arguments)...)};
};


template <class TFunction, class TReturn, class ...TArgs>
concept FunctionWithSignature
= CallableWith<TFunction, TArgs...>
  && requires (TFunction instance, TArgs ... arguments) {
    {instance(std::forward<TArgs>(arguments)...)} -> std::same_as<TReturn>;
};


template <class TFunction, class ...TArgs>
concept AppliableWith = requires (TFunction function, std::tuple<TArgs...> tuple)
{
    {std::apply(function, tuple)};
};


} // namespace cie::concepts


#endif