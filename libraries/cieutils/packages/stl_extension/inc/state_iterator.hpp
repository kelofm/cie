#ifndef CIE_UTILS_STL_EXTENSION_STATE_ITERATOR_HPP
#define CIE_UTILS_STL_EXTENSION_STATE_ITERATOR_HPP

// --- Internal Includes ---
#include "packages/compile_time/packages/concepts/inc/container_concepts.hpp"
#include "packages/types/inc/types.hpp"

// --- STL Includes ---
#include <vector>
#include <iterator>
#include <functional>


namespace cie::utils {


/// @ingroup cieutils
template <class IteratorType>
requires concepts::Container<typename std::iterator_traits<IteratorType>::value_type>
class StateIterator
{
public:
    typedef IteratorType                                                iterator_type;
    typedef typename std::iterator_traits<iterator_type>::value_type    subcontainer_type;
    typedef typename subcontainer_type::const_iterator                  subiterator_type;
    typedef std::vector<subiterator_type>                               state_container;

public:
    StateIterator( const StateIterator& copy ) = delete;
    StateIterator& operator=( const StateIterator& copy ) = delete;

    StateIterator( IteratorType begin,
                   IteratorType end );

    void reset();

    virtual StateIterator& operator++();
    virtual const state_container& operator*() const;

protected:
    explicit StateIterator(); // <-- for initialization of derived classes only

protected:
    IteratorType    _begin;
    IteratorType    _end;
    state_container _state;
}; // class StateIterator





namespace detail {
template <concepts::Container TContainer>
struct container_reference_wrapper : public std::reference_wrapper<TContainer>
{
    using value_type = typename TContainer::value_type;

    using size_type = typename TContainer::size_type;

    using difference_type = typename TContainer::difference_type;

    using iterator = typename TContainer::iterator;

    using const_iterator = typename TContainer::const_iterator;

    container_reference_wrapper(TContainer& container)
        : std::reference_wrapper<TContainer>(container)
    {}

    Size size() const noexcept
    {return this->get().size();}

    iterator begin() noexcept
    {return this->get().begin();}

    const_iterator begin() const noexcept
    {return this->get().begin();}

    iterator end() noexcept
    {return this->get().end();}

    const_iterator end() const noexcept
    {return this->get().end();}

    const_iterator cbegin() const noexcept
    {return this->get().cbegin();}

    const_iterator cend() const noexcept
    {return this->get().cend();}
}; // struct container_reference_wrapper
} // namespace detail



/// @brief Helper class for permutations
/// @ingroup cieutils
template <concepts::Container TContainer>
class InternalStateIterator :
    public StateIterator<typename std::vector<detail::container_reference_wrapper<const TContainer>>::const_iterator>
{
public:
    InternalStateIterator(const TContainer& container,
                          Size size);

protected:
    std::vector<detail::container_reference_wrapper<const TContainer>> _dummyContainer;
};


// ---------------------------------------------------------
// CONVENIENCE FUNCTIONS
// ---------------------------------------------------------

/// @ingroup cieutils
template <class IteratorType>
StateIterator<IteratorType>
makeStateIterator(  IteratorType begin,
                    IteratorType end )
requires concepts::Container<typename std::iterator_traits<IteratorType>::value_type>;


/// @ingroup cieutils
template <concepts::Container ContainerType>
StateIterator<typename ContainerType::const_iterator>
makeStateIterator( const ContainerType& container )
requires concepts::Container<typename ContainerType::value_type>;


/** @details Create a state iterator structure with 'size' components
 *  that has |container.size()|^size number of states.
 *  This basically creates an iterator through all possible
 *  combinations of the container values.
 *  @ingroup cieutils
 */
template <concepts::Container ContainerType>
InternalStateIterator<const ContainerType>
makeInternalStateIterator(const ContainerType& container,
                          Size size);


} // namespace cie::utils

#include "packages/stl_extension/impl/state_iterator_impl.hpp"

#endif
