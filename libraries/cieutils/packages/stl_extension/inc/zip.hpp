#ifndef CIE_UTILS_STL_EXTENSION_ZIP_HPP
#define CIE_UTILS_STL_EXTENSION_ZIP_HPP

// --- Utility Includes ---
#include "packages/types/inc/modifiers.hpp" // CopyConstQualifier
#include "packages/macros/inc/checks.hpp"
#include "packages/compile_time/packages/concepts/inc/container_concepts.hpp"
#include "packages/compile_time/packages/parameter_pack/inc/Size.hpp"

// -- STL Includes ---
#include <tuple>
#include <type_traits>
#include <utility> // std::index_sequence
#include <limits> // std::numeric_limits<T>::max


namespace cie::utils {


namespace detail {


template <class TIterator, class TDistance>
TIterator advanceAndReturn(TIterator it, TDistance distance)
{
    std::advance(it, distance);
    return it;
}


/// @todo Tests fail - fixme.
template <class ...TContainers>
class ZipIterator final
{
public:
    /// Tuple of value references with matching const qualifiers
    using ValueType = std::tuple<
        typename CopyConstQualifier<
            TContainers,
            typename std::remove_reference<TContainers>::type::value_type
        >::Type&...
    >;

    /** Utility class for selecting qualified iterator types
     *  Select TContainer::const_iterator if the input container
     *  is const, or TContainer::iterator otherwise.
     */
    template <class TContainer>
    struct QualifiedIterator
    {};

    template <concepts::Const TContainer>
    struct QualifiedIterator<TContainer>
    {using Type = typename std::remove_reference<TContainer>::type::const_iterator;};

    template <concepts::NonConst TContainer>
    struct QualifiedIterator<TContainer>
    {using Type = typename std::remove_reference<TContainer>::type::iterator;};

    /// Tuple of const or non-const iterators depending on container qualifiers
    using IteratorTuple = std::tuple<typename QualifiedIterator<TContainers>::Type...>;

public:
    ZipIterator(IteratorTuple&& r_it)
        : _it(r_it)
    {}

    ZipIterator(ZipIterator&& r_rhs) = default;

    ZipIterator(const ZipIterator& r_rhs) = default;

    ZipIterator& operator=(ZipIterator&& r_rhs) = default;

    ZipIterator& operator=(const ZipIterator& r_rhs) = default;

    ValueType operator*()
    {return this->makeRefTuple<ValueType>();}

    const ValueType operator*() const
    {return this->makeRefTuple<const ValueType>();}

    ZipIterator& operator++()
    {
        IteratorManipulator<std::tuple_size<IteratorTuple>::value-1>::increment(_it);
        return *this;
    }

    ZipIterator operator++(int)
    {
        ZipIterator copy(*this);
        ++(*this);
        return copy;
    }

    ZipIterator& operator--()
    {
        IteratorManipulator<std::tuple_size<IteratorTuple>::value-1>::decrement(_it);
        return *this;
    }

    ZipIterator operator--(int)
    {
        ZipIterator copy(*this);
        --(*this);
        return copy;
    }

    ZipIterator& operator+=(int distance)
    {
        IteratorManipulator<sizeof...(TContainers)-1>::advance(_it, distance);
        return *this;
    }

    ZipIterator& operator-=(int distance)
    {
        IteratorManipulator<sizeof...(TContainers)-1>::advance(_it, distance);
        return *this;
    }

    bool operator==(const ZipIterator& r_rhs) const
    {return IteratorManipulator<sizeof...(TContainers)-1>::isEqual(_it, r_rhs._it);}

    bool operator<(const ZipIterator& r_rhs) const
    {return IteratorManipulator<sizeof...(TContainers)-1>::isLess(_it, r_rhs._it);}

private:
    template <Size Index>
    struct IteratorManipulator
    {
        static void advance(IteratorTuple& r_it, int distance)
        {
            std::advance(std::get<Index>(r_it), distance);
            if constexpr (0 < Index) IteratorManipulator<Index-1>::advance(r_it, distance);
        };

        static void increment(IteratorTuple& r_it)
        {
            ++std::get<Index>(r_it);
            if constexpr (0 < Index) IteratorManipulator<Index-1>::increment(r_it);
        }

        static void decrement(IteratorTuple& r_it)
        {
            --std::get<Index>(r_it);
            if constexpr (0 < Index) IteratorManipulator<Index-1>::increment(r_it);
        }

        static bool isEqual(const IteratorTuple& r_lhs, const IteratorTuple& r_rhs)
        {
            // Assume that one check yields the same result for all stored iterators
            return std::get<Index>(r_lhs) == std::get<Index>(r_rhs);
        }

        static bool isLess(const IteratorTuple& r_lhs, const IteratorTuple& r_rhs)
        {
            // Assume that one check yields the same result for all stored iterators
            return std::get<Index>(r_lhs) < std::get<Index>(r_rhs);
        }
    }; // struct IteratorManipulator

    template <class TValue>
    TValue makeRefTuple();

    template <class TValue, std::size_t ...Indices>
    TValue makeRefTupleImpl(std::index_sequence<Indices...> sequence);

private:
    IteratorTuple _it;
}; // class ZipIterator


template <class ...TContainers>
class ZipObject final
{
public:
    using Iterator = ZipIterator<TContainers...>;

    using ValueType = typename Iterator::ValueType;

    using iterator = Iterator;

    using value_type = ValueType;

public:
    ZipObject(TContainers&&... r_containers)
        : _containers(r_containers...),
          _size(MinimumSize<sizeof...(TContainers)-1>::get(_containers))
    {}

    ZipObject(Size size, TContainers&&... r_containers)
        : _containers(r_containers...),
          _size(size)
    {}

    Iterator begin()
    {return this->beginImpl(std::make_index_sequence<sizeof...(TContainers)>());}

    Iterator end()
    {return this->endImpl(std::make_index_sequence<sizeof...(TContainers)>());}

private:
    using ContainerTuple = std::tuple<TContainers...>;

    template <Size Index>
    struct MinimumSize
    {
        static Size get(const ContainerTuple& r_containers)
        {return MinimumSize<Index>::getImpl(std::numeric_limits<Size>::max(), r_containers);}

        static Size getImpl(Size min, const ContainerTuple& r_containers)
        {
            const Size current = std::get<Index>(r_containers).size();
            const Size newMin = min < current ? min : current;
            if constexpr (0 < Index) return MinimumSize<Index-1>::getImpl(newMin, r_containers);
            else return newMin;
        }
    }; // struct MinimumSize

    template <std::size_t ...Indices>
    Iterator beginImpl(std::index_sequence<Indices...> sequence);

    template <std::size_t ...Indices>
    Iterator endImpl(std::index_sequence<Indices...> sequence);

    ContainerTuple _containers;

    const Size _size;
}; // class ZipObject


template <class ...Ts>
struct ZippableHelper
{
    using Tuple = std::tuple<Ts...>;

    template <Size I>
    static constexpr bool isZippable()
    {
        bool output = concepts::Container<typename std::decay<typename std::tuple_element<I,Tuple>::type>::type>;
        if constexpr (I == 0)
            return output;
        else
            return output && ZippableHelper<Ts...>::isZippable<I-1>();
    }
}; // struct ZippableHelper


template <class ...Ts>
concept Zippables
= ZippableHelper<Ts...>::template isZippable<ct::PackSize<Ts...>-1>();


} // namespace detail


/// @ingroup cieutils
template <detail::Zippables ...TContainers>
detail::ZipObject<TContainers...> zip(Size size, TContainers&&... r_containers)
{
    return detail::ZipObject<TContainers...>(size, std::forward<TContainers>(r_containers)...);
}


/// @ingroup cieutils
template <detail::Zippables ...TContainers>
detail::ZipObject<TContainers...> zip(TContainers&&... r_containers)
{
    return detail::ZipObject<TContainers...>(std::forward<TContainers>(r_containers)...);
}


} // namespace cie::utils

#include "packages/stl_extension/impl/zip_impl.hpp"

#endif
