#pragma once

// --- Internal Includes ---
#include "packages/compile_time/packages/concepts/inc/iterator_concepts.hpp"


namespace cie::concepts {


namespace detail {

// Member requirements
template < class T,
           class TArgument = typename T::size_type,
           class TReturn = typename T::reference >
concept HasAt
= requires (T instance, TArgument argument)
{
    {instance.at(argument)} -> std::same_as<TReturn>;
};


template < class T,
           class TArgument = typename T::size_type,
           class TReturn = typename T::reference >
concept HasAccessOperator
= requires (T instance, TArgument argument)
{
    {instance[argument]} -> std::same_as<TReturn>;
};


template < class T,
           class TReturn = typename T::reference >
concept HasFront
= requires (T instance)
{
    {instance.front()} -> std::same_as<TReturn>;
};


template < class T,
           class TReturn = typename T::reference >
concept HasBack
= requires (T instance)
{
    {instance.back()} -> std::same_as<TReturn>;
};


template < class T,
           class ...TArgs>
concept HasResize
= requires (T instance, TArgs... arguments)
{
    {instance.resize(arguments...)} -> std::same_as<void>;
};


template < class T,
           class TArgument = typename T::size_type,
           class TReturn = void >
concept HasReserve
= requires (T instance, TArgument argument)
{
    {instance.reserve(argument)} -> std::same_as<TReturn>;
};


template <class T>
concept HasClear
= requires (T instance)
{
    {instance.clear()};
};


template < class T,
           class TReturn = typename T::iterator >
concept HasErase
= requires (T instance)
{
    {instance.erase(instance.begin(), instance.end())} -> std::same_as<TReturn>;
};


template < class T,
           class TArgument = const typename T::value_type&,
           class TReturn = void >
concept HasPushFront
= requires (T instance, TArgument argument)
{
    {instance.push_front(argument)} -> std::same_as<TReturn>;
};


template < class T,
           class TArgument = const typename T::value_type&,
           class TReturn = void >
concept HasPushBack
= requires (T instance, TArgument argument)
{
    {instance.push_back(argument)} -> std::same_as<TReturn>;
};


template < class T,
           class ...TArgs >
concept HasEmplaceFront
= requires (T instance, TArgs&&... rArgs)
{
    {instance.emplace_front(std::forward<TArgs>(rArgs)...)} -> std::same_as<typename T::reference>;
};


template < class T,
           class ...TArgs >
concept HasEmplaceBack
= requires (T instance, TArgs&&... rArgs)
{
    {instance.emplace_back(std::forward<TArgs>(rArgs)...)} -> std::same_as<typename T::reference>;
};


template <class T, class ...TArgs>
concept HasInsert
= requires (T instance, TArgs&&... rTArgs) {
    {instance.insert(std::forward<TArgs>(rTArgs)...)};
};


template < class T,
           class TReturn = void >
concept HasPopFront
= requires (T instance)
{
    {instance.pop_back()} -> std::same_as<TReturn>;
};


template < class T,
           class TReturn = void >
concept HasPopBack
= requires (T instance)
{
    {instance.pop_back()} -> std::same_as<TReturn>;
};


template < class T,
           class TReturn = void >
concept HasSwap
= requires (T instance, T swap)
{
    {instance.swap(swap)} -> std::same_as<TReturn>;
};


template <class T>
concept HasData
= requires (T instance, const T constInstance)
{
    {instance.data()} -> std::same_as<typename T::value_type*>;
    {constInstance.data()} -> std::same_as<const typename T::value_type*>;
};


} // namespace detail


template <class T, class TValue = void>
concept Container
= requires (T instance, const T constInstance)
{
    typename T::value_type;
    typename T::size_type;
    typename T::iterator;
    //typename T::const_iterator;
    {instance.size()};
    {instance.begin()};
    {instance.end()};
    {constInstance.begin()};
    {constInstance.end()};
} && (std::is_same_v<TValue,void> || std::is_same_v<typename T::value_type,TValue>);


template <class TContainer>
concept NumericContainer
= Container<typename std::decay<TContainer>::type>
  && Numeric<typename std::decay<TContainer>::type::value_type>;


template <class TContainer, class TValue = void>
concept PointerContainer
= Container<typename std::remove_reference<TContainer>::type>
  && Pointer<typename std::remove_reference<TContainer>::type::value_type, TValue>;


template <class TContainer>
concept NonPointerContainer
= Container<typename std::remove_reference<TContainer>::type>
  && !PointerContainer<TContainer>;


template <class TContainer>
concept IteratorContainer
= Container<typename std::remove_reference<TContainer>::type>
  && Iterator<typename std::remove_reference<TContainer>::type::value_type>;


template <class TContainer, class InterfaceType>
concept InterfaceContainer
= Container<typename std::remove_reference<TContainer>::type>
  && Pointer<typename std::remove_reference<TContainer>::type::value_type>
  && DerivedFrom<typename std::pointer_traits<typename std::remove_reference<TContainer>::type::value_type>::element_type,InterfaceType>;


// ---------------------------------------------------------
// SPECIALIZED STL CONTAINERS
// ---------------------------------------------------------

template <class TContainer>
concept ResizableContainer
= Container<TContainer> && detail::HasResize<TContainer,std::size_t>;


template <class TContainer>
concept ReservableContainer
= Container<TContainer>
  && requires (TContainer instance)
{
    {instance.reserve(1)};
};


// TODO: improve definition
template <class TContainer>
concept StaticContainer
=   Container<TContainer>
    && !ResizableContainer<TContainer>;


} // namespace cie::concepts
