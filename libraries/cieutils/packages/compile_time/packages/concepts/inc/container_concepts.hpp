#ifndef CIE_CIEUTILS_CONTAINER_CONCEPTS_HPP
#define CIE_CIEUTILS_CONTAINER_CONCEPTS_HPP

// --- Internal Includes ---
#include "packages/compile_time/packages/concepts/inc/iterator_concepts.hpp"


namespace cie::concepts {


namespace detail {

// Member requirements
template < class T,
           class ArgumentType = typename T::size_type,
           class ReturnType = typename T::reference >
concept HasAt
= requires (T instance, ArgumentType argument)
{
    {instance.at(argument)} -> std::same_as<ReturnType>;
};


template < class T,
           class ArgumentType = typename T::size_type,
           class ReturnType = typename T::reference >
concept HasAccessOperator
= requires (T instance, ArgumentType argument)
{
    {instance[argument]} -> std::same_as<ReturnType>;
};


template < class T,
           class ReturnType = typename T::reference >
concept HasFront
= requires (T instance)
{
    {instance.front()} -> std::same_as<ReturnType>;
};


template < class T,
           class ReturnType = typename T::reference >
concept HasBack
= requires (T instance)
{
    {instance.back()} -> std::same_as<ReturnType>;
};


template < class T,
           class ...Arguments>
concept HasResize
= requires (T instance, Arguments... arguments)
{
    {instance.resize(arguments...)} -> std::same_as<void>;
};


template < class T,
           class ArgumentType = typename T::size_type,
           class ReturnType = void >
concept HasReserve
= requires (T instance, ArgumentType argument)
{
    {instance.reserve(argument)} -> std::same_as<ReturnType>;
};


template <class T>
concept HasClear
= requires (T instance)
{
    {instance.clear()};
};


template < class T,
           class ReturnType = typename T::iterator >
concept HasErase
= requires (T instance)
{
    {instance.erase(instance.begin(), instance.end())} -> std::same_as<ReturnType>;
};


template < class T,
           class ArgumentType = const typename T::value_type&,
           class ReturnType = void >
concept HasPushFront
= requires (T instance, ArgumentType argument)
{
    {instance.push_front(argument)} -> std::same_as<ReturnType>;
};


template < class T,
           class ArgumentType = const typename T::value_type&,
           class ReturnType = void >
concept HasPushBack
= requires (T instance, ArgumentType argument)
{
    {instance.push_back(argument)} -> std::same_as<ReturnType>;
};


template < class T,
           class ...Arguments >
concept HasEmplaceFront
= requires (T instance, Arguments&&... r_arguments)
{
    {instance.emplace_front(std::forward<Arguments>(r_arguments)...)} -> std::same_as<typename T::reference>;
};


template < class T,
           class ...Arguments >
concept HasEmplaceBack
= requires (T instance, Arguments&&... r_arguments)
{
    {instance.emplace_back(std::forward<Arguments>(r_arguments)...)} -> std::same_as<typename T::reference>;
};


template < class T,
           class ReturnType = void >
concept HasPopFront
= requires (T instance)
{
    {instance.pop_back()} -> std::same_as<ReturnType>;
};


template < class T,
           class ReturnType = void >
concept HasPopBack
= requires (T instance)
{
    {instance.pop_back()} -> std::same_as<ReturnType>;
};


template < class T,
           class ReturnType = void >
concept HasSwap
= requires (T instance, T swap)
{
    {instance.swap(swap)} -> std::same_as<ReturnType>;
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
    typename T::const_iterator;
    {instance.size()};
    {instance.begin()};
    {instance.end()};
    {constInstance.begin()};
    {constInstance.end()};
} && (std::is_same_v<TValue,void> || std::is_same_v<typename T::value_type,TValue>);


template <class ContainerType>
concept NumericContainer
= Container<typename std::decay<ContainerType>::type>
  && Numeric<typename std::decay<ContainerType>::type::value_type>;


template <class TContainer, class TValue = void>
concept PointerContainer
= Container<typename std::remove_reference<TContainer>::type>
  && Pointer<typename std::remove_reference<TContainer>::type::value_type, TValue>;


template <class TContainer>
concept NonPointerContainer
= Container<typename std::remove_reference<TContainer>::type>
  && !PointerContainer<TContainer>;


template <class ContainerType>
concept IteratorContainer
= Container<typename std::remove_reference<ContainerType>::type>
  && Iterator<typename std::remove_reference<ContainerType>::type::value_type>;


template <class ContainerType, class InterfaceType>
concept InterfaceContainer
= Container<typename std::remove_reference<ContainerType>::type>
  && Pointer<typename std::remove_reference<ContainerType>::type::value_type>
  && DerivedFrom<typename std::pointer_traits<typename std::remove_reference<ContainerType>::type::value_type>::element_type,InterfaceType>;


// ---------------------------------------------------------
// SPECIALIZED STL CONTAINERS
// ---------------------------------------------------------

template <class ContainerType>
concept ResizableContainer
= Container<ContainerType>
  && requires (ContainerType instance)
{
    {instance.resize(1)};
};


template <class ContainerType>
concept ReservableContainer
= Container<ContainerType>
  && requires (ContainerType instance)
{
    {instance.reserve(1)};
};


// TODO: improve definition
template <class ContainerType>
concept StaticContainer
=   Container<ContainerType>
    && !ResizableContainer<ContainerType>;


} // namespace cie::concepts

#endif