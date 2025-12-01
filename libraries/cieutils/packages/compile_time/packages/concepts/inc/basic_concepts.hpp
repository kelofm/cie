#ifndef CIE_CIEUTILS_BASIC_CONCEPTS_HPP
#define CIE_CIEUTILS_BASIC_CONCEPTS_HPP

// --- Utility Includes ---
#include "packages/types/inc/types.hpp"

// --- STL Includes ---
#include <concepts>
#include <type_traits>
#include <cstdint>
#include <utility>


namespace cie::concepts {


// ---------------------------------------------------------
// STL Extension
// ---------------------------------------------------------

template <class T, class TT>
concept same_as_noqualifiers
= std::is_same_v<T,TT> || std::is_same_v<typename std::remove_const<T>::type,TT>;

template <class T>
concept Const
= std::is_const_v<typename std::remove_reference<T>::type>;


template <class T>
concept NonConst
= !Const<T>;

template <class T, class ...TT>
concept AnyOf
= (std::is_same_v<T,TT> || ...);

template <class T, class ...TT>
concept NoneOf
= !AnyOf<T,TT...>;

namespace detail {
template <class T> struct AnyNoneHelper {};

template <class ...Ts> struct AnyNoneHelper<std::tuple<Ts...>>
{
    template <class T>
    static constexpr inline bool AnyOfValue = AnyOf<T,Ts...>;

    template <class T>
    static constexpr inline bool NoneOfValue = NoneOf<T,Ts...>;
};
} // namespace detail

template <class T, class TTuple>
concept AnyOfTuple
= detail::AnyNoneHelper<TTuple>::template AnyOfValue<T>;

template <class T, class TTuple>
concept NoneOfTuple
= detail::AnyNoneHelper<TTuple>::template NoneOfValue<T>;

// ---------------------------------------------------------
// Life Cycle
// ---------------------------------------------------------

template <class T>
concept Constructible
= std::is_constructible_v<T>;

template <class T>
concept DefaultConstructible
= std::is_default_constructible_v<T>;

template <class T>
concept MoveConstructible
= std::is_move_constructible_v<T>;

template <class T>
concept CopyConstructible
= std::is_copy_constructible_v<T>;

template <class T>
concept MoveAssignable
= std::is_move_assignable_v<T>;

template <class T>
concept CopyAssignable
= std::is_copy_assignable_v<T>;

template <class T>
concept Destructible
= std::is_destructible_v<T>;

template <class T, class ...Arguments>
concept ConstructibleFrom
= requires(Arguments&&... r_arguments)
{
    {T(std::forward<Arguments>(r_arguments)...)};
};


// ---------------------------------------------------------
// Polymorphism
// ---------------------------------------------------------

namespace detail {
template<typename T>
struct void_ {typedef void type;};

template<typename T, typename = void>
struct CanBeDerivedFrom {
  static const bool value = false;
};

template<typename T>
struct CanBeDerivedFrom<T, typename void_<int T::*>::type> {
  static const bool value = true;
};
} // namespace detail


template <class T>
concept Deriveable
= detail::CanBeDerivedFrom<T>::value
    && !std::is_integral_v<T>
    && !std::is_floating_point_v<T>;


template <class T>
concept NonDeriveable
= !Deriveable<T>;

template <class DerivedType, class BaseType>
concept DerivedFrom
= std::derived_from<DerivedType,BaseType>;


// ---------------------------------------------------------
// Numerics
// ---------------------------------------------------------

template <class T>
concept Integral
= std::is_integral_v<T>;

template <class T>
concept Incrementable
= requires (T instance)
{
    {++instance};
    {instance++};
};


template <class T>
concept Decrementable
= requires (T instance)
{
    {--instance};
    {instance--};
};


template <class T>
concept Addable
= requires (T instance)
{
    {instance + instance}     -> std::same_as<T>;
    {instance += instance}    -> std::same_as<T&>;
};


template <class TLeft, class TRight = TLeft>
concept Subtractable
= requires (TLeft left, TRight right)
{
    {left - right}     -> std::same_as<TLeft>;
    {left -= right}    -> std::same_as<TLeft&>;
};


template <class T>
concept Multiplicable
= requires (T instance)
{
    {instance * instance}     -> std::same_as<T>;
    {instance *= instance}    -> std::same_as<T&>;
};


template <class T>
concept Divisible
= requires (T instance)
{
    {instance / instance}     -> std::same_as<T>;
    {instance /= instance}    -> std::same_as<T&>;
};


template <class T>
concept Numeric
=  Incrementable<T>
&& Decrementable<T>
&& Addable<T>
&& Subtractable<T>
&& Multiplicable<T>
&& Divisible<T>;


template <class T>
concept UnsignedInteger
= AnyOf<T, uint8_t,
           uint16_t,
           uint32_t,
           uint64_t,
           uint_fast8_t,
           uint_fast16_t,
           uint_fast32_t,
           uint_fast64_t,
           uint_least8_t,
           uint_least16_t,
           uint_least32_t,
           uint_least64_t,
           unsigned short,
           unsigned int,
           unsigned long,
           unsigned long long>;

template <class T>
concept SignedInteger
= AnyOf<T, int8_t,
           int16_t,
           int32_t,
           int64_t,
           int_fast8_t,
           int_fast16_t,
           int_fast32_t,
           int_fast64_t,
           int_least8_t,
           int_least16_t,
           int_least32_t,
           int_least64_t>;

template <class T>
concept Integer
= UnsignedInteger<T> || SignedInteger<T>;


// ---------------------------------------------------------
// Casting and Conversion
// ---------------------------------------------------------

template <class SourceType, class TargetType>
concept ConvertibleTo
= requires (SourceType instance)
{
    {TargetType(instance)};
};


// ---------------------------------------------------------
// Pointers
// ---------------------------------------------------------

template <class TDereferencable, class TDereferenced = void>
concept Dereferencable
= (std::is_same_v<TDereferenced,void> && requires (TDereferencable instance)
{
    {*instance};
    {instance.operator->()};
}) || requires (TDereferencable instance)
{
    {*instance}                 -> std::same_as<Ref<TDereferenced>>;
    {instance.operator->()}     -> std::same_as<Ptr<TDereferenced>>;
};


template <class TPointer, class TElement = void>
concept RawPointer
= std::is_pointer_v<TPointer> && (std::is_same_v<TElement,void> || std::is_same_v<typename std::pointer_traits<TPointer>::element_type,TElement>);


template <class TPointer, class TElement = void>
concept Pointer
= RawPointer<TPointer,TElement> || Dereferencable<TPointer,TElement>;


template <class TPointer, class TElement>
concept NonRawPointer
= Pointer<TPointer,TElement> && !RawPointer<TPointer,TElement>;


template <class T>
concept NonPointer
= !Pointer<T>;


// ---------------------------------------------------------
// Memory
// ---------------------------------------------------------


template <class TAllocator, class TValue>
concept Allocator
= requires (TAllocator instance,
            std::size_t count,
            TValue* pValue) {
    //typename TAllocator::pointer;
    //typename TAllocator::const_pointer;
    //typename TAllocator::void_pointer;
    //typename TAllocator::const_void_pointer;
    typename TAllocator::value_type;
    //typename TAllocator::size_type;
    //typename TAllocator::difference_type;
    {instance.allocate(count)} -> std::same_as<TValue*>;
    {instance.deallocate(pValue, count)};
    //{instance.max_size()} -> std::same_as<typename TAllocator::size_type>;
};


// ---------------------------------------------------------
// MISC
// ---------------------------------------------------------

namespace detail {
template <class T>
struct FalseType : std::false_type {};
} // namespace detail


template <class T>
concept Unsatisfiable
= detail::FalseType<T>::value;

} // namespace cie::concepts

#endif