#ifndef CIE_UTILS_UNSTRUCTURED_VIEW_HPP
#define CIE_UTILS_UNSTRUCTURED_VIEW_HPP

// --- Utility Includes ---
#include "packages/stl_extension/inc/DynamicArray.hpp"
#include "packages/compile_time/packages/concepts/inc/container_concepts.hpp"
#include "packages/macros/inc/typedefs.hpp"

// --- STL Includes ---
#include <iterator>
#include <memory>


namespace cie::utils {


///@addtogroup cieutils
///@{

/// @brief A container that stores pointers to items but provides access through references.
template <class TValue, class TContainer = DynamicArray<TValue*>>
requires (concepts::Container<TContainer,const TValue*> || concepts::Container<TContainer,TValue*>)
class UnstructuredView
{
private:
    using Container = TContainer;

    inline constexpr static const bool isMutableView = !std::is_const_v<typename std::pointer_traits<typename TContainer::value_type>::element_type>;

    template <class T, class TWrapped>
    class Iterator
    {
    public:
        using value_type = T;

        using pointer = value_type*;

        using const_pointer = const value_type*;

        using reference = value_type&;

        using const_reference = const value_type&;

        using difference_type = std::ptrdiff_t;

        using iterator_category = std::random_access_iterator_tag;

    private:
        friend class UnstructuredView;

    public:
        Iterator() = delete;

        CIE_DEFINE_CLASS_DEFAULT_MOVES(Iterator)

        CIE_DEFINE_CLASS_DEFAULT_COPIES(Iterator)

        reference operator*() requires isMutableView {return **_it;}

        const_reference operator*() const {return **_it;};

        pointer operator->() requires isMutableView {return *_it;}

        const_pointer operator->() const {return *_it;}

        Iterator& operator++() {++_it; return *this;}

        Iterator operator++(int) {Iterator copy(*this); ++_it; return copy;}

        Iterator& operator--() {--_it; return *this;}

        Iterator operator--(int) {Iterator copy(*this); --_it; return copy;}

        Iterator& operator+=(difference_type difference) {_it+=difference; return *this;}

        Iterator& operator-=(difference_type difference) {_it-=difference; return *this;}

        friend Iterator operator+(Iterator lhs, difference_type rhs) {lhs._it+=rhs; return lhs;}

        friend Iterator operator+(difference_type lhs, Iterator rhs) {rhs._it+=lhs; return rhs;}

        friend Iterator operator-(Iterator lhs, difference_type rhs) {lhs._it-=rhs; return lhs;}

        friend difference_type operator-(Iterator lhs, Iterator rhs) {return lhs._it - rhs._it;}

        friend bool operator<(Iterator lhs, Iterator rhs) {return lhs._it < rhs._it;}

        friend bool operator!=(Iterator lhs, Iterator rhs) {return lhs._it != rhs._it;}

    private:
        Iterator(TWrapped it) : _it(it) {}

    private:
        TWrapped _it;
    }; // class Iterator

    using ContainerValue = typename std::pointer_traits<typename TContainer::value_type>::element_type;

public:
    using value_type = TValue;

    using reference = TValue&;

    using const_reference = const TValue&;

    using iterator = Iterator<ContainerValue,typename TContainer::iterator>;

    using const_iterator = Iterator<const ContainerValue,typename TContainer::iterator>;

    using reverse_iterator = Iterator<ContainerValue,typename TContainer::reverse_iterator>;

    using const_reverse_iterator = Iterator<const ContainerValue,typename TContainer::reverse_iterator>;

    using size_type = Size;

    CIE_DEFINE_CLASS_POINTERS(UnstructuredView)

public:
    CIE_DEFINE_CLASS_DEFAULTS(UnstructuredView)

    /// @brief Construct directly from a moved container of pointers.
    UnstructuredView(TContainer&& r_pointers);

    /// @brief Construct from a range of pointers.
    template <class TIterator>
    requires (concepts::Iterator<TIterator,TValue*> || concepts::Iterator<TIterator,const TValue*>)
    UnstructuredView(TIterator begin, TIterator end);

    /// @brief Construct from a range of instances by copying their addresses.
    template <class TIterator>
    requires (concepts::Iterator<TIterator,const TValue> || concepts::Iterator<TIterator,TValue>)
    UnstructuredView(TIterator begin, TIterator end);

    const_reference at(size_type index) const;

    reference at(size_type index)
    requires isMutableView;

    const_reference operator[](size_type index) const;

    reference operator[](size_type index)
    requires isMutableView;

    void push_back(ContainerValue& r_value);

    void emplace_back(ContainerValue& r_value);

    void reserve(size_type capacity);

    /// @brief Get the number of stored pointers.
    Size size() const noexcept;

    Size capacity() const noexcept;

    const_iterator begin() const noexcept;

    iterator begin() noexcept
    requires isMutableView;

    const_iterator end() const noexcept;

    iterator end() noexcept
    requires isMutableView;

    const_reverse_iterator rbegin() const noexcept;

    reverse_iterator rbegin() noexcept
    requires isMutableView;

    const_reverse_iterator rend() const noexcept;

    reverse_iterator rend() noexcept
    requires isMutableView;

private:
    mutable TContainer _pointers;
}; // class UnstructuredView

///@}

} // namespace cie::utils

#include "packages/ranges/impl/UnstructuredView_impl.hpp"

#endif
