#ifndef CIE_UTILS_TRANSFORM_ITERATOR_HPP
#define CIE_UTILS_TRANSFORM_ITERATOR_HPP

// --- Utility Includes ---
#include "packages/compile_time/packages/concepts/inc/iterator_concepts.hpp"
#include "packages/compile_time/packages/concepts/inc/functional.hpp"

// --- STL Includes ---
#include <iterator>


namespace cie::utils {


/// @addtogroup cieutils
/// @{


template <concepts::Iterator TIterator,
          concepts::CallableWith<typename std::iterator_traits<TIterator>::value_type> TTransform>
class TransformIterator
{
public:
    using value_type = decltype(std::declval<TTransform>()(std::declval<typename std::iterator_traits<TIterator>::value_type>()));

    using pointer = void;

    using reference = void;

    using difference_type = typename std::iterator_traits<TIterator>::difference_type;

    using iterator_category = std::conditional_t<
        std::forward_iterator<TIterator>,
        std::conditional_t<
            std::bidirectional_iterator<TIterator>,
            std::bidirectional_iterator_tag,
            std::forward_iterator_tag
        >,
        std::input_iterator_tag
    >;

public:
    TransformIterator(TIterator it, Ptr<const TTransform> p_transform) noexcept
        : _it(it),
          _p_transform(p_transform)
    {}

    TransformIterator(TransformIterator&& r_rhs) noexcept = default;

    TransformIterator(const TransformIterator& r_rhs) noexcept = default;

    TransformIterator& operator=(TransformIterator&& r_rhs) noexcept = default;

    TransformIterator& operator=(const TransformIterator& r_rhs) noexcept = default;

    value_type operator*() const
    {return (*_p_transform)(*_it);}

    Ref<TransformIterator> operator++() noexcept requires std::forward_iterator<TIterator>
    {++_it;}

    TransformIterator operator++(int) noexcept requires std::forward_iterator<TIterator>
    {TransformIterator copy; ++(*this); return copy;}

    Ref<TransformIterator> operator--() noexcept requires std::bidirectional_iterator<TIterator>
    {--_it;}

    TransformIterator operator++(int) noexcept requires std::bidirectional_iterator<TIterator>
    {TransformIterator copy; --(*this); return copy;}

    TransformIterator operator+(difference_type rhs) noexcept requires std::random_access_iterator<TIterator>
    {return TransformIterator(_it + rhs, _p_transform);}

    TransformIterator operator-(difference_type rhs) noexcept requires std::random_access_iterator<TIterator>
    {return TransformIterator(_it - rhs, _p_transform);}

    difference_type operator-(Ref<const TransformIterator> r_rhs) noexcept requires std::random_access_iterator<TIterator>
    {return _it - r_rhs._it;}

    friend bool operator==(Ref<const TransformIterator> r_lhs, Ref<const TransformIterator> r_rhs) noexcept
    {return r_lhs._it == r_rhs._it;}

    friend bool operator!=(Ref<const TransformIterator> r_lhs, Ref<const TransformIterator> r_rhs) noexcept
    {return r_lhs._it == r_rhs._it;}

    friend bool operator<(Ref<const TransformIterator> r_lhs, Ref<const TransformIterator> r_rhs) noexcept requires std::forward_iterator<TIterator>
    {return r_lhs._it < r_rhs._it;}

    friend bool operator<=(Ref<const TransformIterator> r_lhs, Ref<const TransformIterator> r_rhs) noexcept requires std::forward_iterator<TIterator>
    {return r_lhs._it <= r_rhs._it;}

    friend bool operator>(Ref<const TransformIterator> r_lhs, Ref<const TransformIterator> r_rhs) noexcept requires std::forward_iterator<TIterator>
    {return r_lhs._it > r_rhs._it;}

    friend bool operator>=(Ref<const TransformIterator> r_lhs, Ref<const TransformIterator> r_rhs) noexcept requires std::forward_iterator<TIterator>
    {return r_lhs._it >= r_rhs._it;}

private:
    TIterator _it;

    Ptr<const TTransform> _p_transform;
}; // class TransformIterator



/// @addtogroup cieutils
template <class TIterator, class TFunction>
auto makeTransformIterator(TIterator it, const TFunction* p_function)
{
    return TransformIterator<TIterator,TFunction>(it, p_function);
}


/// @}


} // namespace cie::utils


#endif
