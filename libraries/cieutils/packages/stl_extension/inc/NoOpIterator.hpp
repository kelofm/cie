#ifndef CIE_UTILS_NOOP_ITERATOR_HPP
#define CIE_UTILS_NOOP_ITERATOR_HPP

// --- Utility Includes ---
#include "packages/compile_time/packages/concepts/inc/basic_concepts.hpp"
#include "packages/compile_time/packages/concepts/inc/iterator_concepts.hpp"
#include "packages/types/inc/types.hpp"


namespace cie::utils {


template <concepts::Iterator T>
class NoOpIterator
{
public:
    using value_type = T;

    using pointer = Ptr<T>;

    using reference = Ref<T>;

    using difference_type = typename T::difference_type;

    using iterator_category = typename T::iterator_category;

public:
    NoOpIterator() noexcept = default;

    NoOpIterator(T it) noexcept : _it(it) {}

    reference operator*() const noexcept {return _it;}

    pointer operator->() const noexcept {return &_it;}

    NoOpIterator& operator++() noexcept {++_it; return *this;}

    NoOpIterator operator++(int) noexcept {NoOpIterator copy=*this; ++(*this); return copy;}

    NoOpIterator& operator--() noexcept {--_it; return *this;}

    NoOpIterator operator--(int) noexcept {NoOpIterator copy=*this; --(*this); return copy;}

    NoOpIterator& operator+=(difference_type difference) noexcept {_it+=difference; return *this;}

    NoOpIterator& operator-=(difference_type difference) noexcept
    requires concepts::Subtractable<T,difference_type>
    {_it-=difference; return *this;}

    NoOpIterator operator+(difference_type difference) const noexcept {return _it+difference;}

    NoOpIterator operator-(difference_type difference) const noexcept
    requires concepts::Subtractable<T,difference_type>
    {return _it-difference;}

    difference_type operator-(NoOpIterator rhs) const noexcept
    requires concepts::Subtractable<T,T>
    {return _it - rhs._it;}

    bool operator==(NoOpIterator rhs) const noexcept {return _it==rhs._it;}

    bool operator!=(NoOpIterator rhs) const noexcept {return _it!=rhs._it;}

    bool operator<(NoOpIterator rhs) const noexcept {return _it<rhs._it;}

    bool operator>(NoOpIterator rhs) const noexcept {return _it>rhs._it;}

    bool operator<=(NoOpIterator rhs) const noexcept {return _it<=rhs._it;}

    bool operator>=(NoOpIterator rhs) const noexcept {return _it>=rhs._it;}

private:
    mutable T _it;
};


template <class T>
NoOpIterator<T> makeNoOpIterator(T it) {return it;}


} // namespace cie::utils


#endif
