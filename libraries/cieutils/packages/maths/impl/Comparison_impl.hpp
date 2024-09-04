#ifndef CIE_UTILS_COMPARISON_IMPL_HPP
#define CIE_UTILS_COMPARISON_IMPL_HPP

// --- Utility Includes ---
#include "packages/maths/inc/Comparison.hpp"

// --- STL Includes ---
#include <limits>
#include <cmath>
#include <algorithm> // all_of


namespace cie::utils {


template <class TInteger>
bool IntegerComparison<TInteger>::equal(const TInteger left,
                                        const TInteger right) const noexcept
{
    return left == right;
}


template <class TInteger>
bool IntegerComparison<TInteger>::less(const TInteger left,
                                       const TInteger right) const noexcept
{
    return left < right;
}


template <class TFloat>
FloatComparison<TFloat>::FloatComparison() noexcept
    : FloatComparison(std::numeric_limits<TFloat>::min(),
                      std::numeric_limits<TFloat>::epsilon())
{
}


template <class TFloat>
FloatComparison<TFloat>::FloatComparison(const TFloat absoluteTolerance,
                                         const TFloat relativeTolerance) noexcept
    : _absoluteTolerance(absoluteTolerance),
      _relativeTolerance(relativeTolerance)
{
}


template <class TFloat>
bool FloatComparison<TFloat>::equal(const TFloat left,
                                    const TFloat right) const noexcept
{
    const TFloat norm = std::min(std::abs(left) + std::abs(right),
                                 std::numeric_limits<TFloat>::max());
    return std::abs(left - right) < std::max(_absoluteTolerance, _relativeTolerance * norm);
}


template <class TFloat>
bool FloatComparison<TFloat>::less(const TFloat left,
                                   const TFloat right) const noexcept
{
    return left < right && !this->equal(left, right);
}


template <class TContainer>
requires concepts::Container<TContainer>
template <class ...TArgs>
Comparison<TContainer>::Comparison(TArgs&&... rArgs) noexcept
    : _componentComparison(std::forward<TArgs>(rArgs)...)
{
}


template <class TContainer>
requires concepts::Container<TContainer>
bool Comparison<TContainer>::equal(Ref<const TContainer> rLeft,
                                   Ref<const TContainer> rRight) const noexcept
{
    bool output = rLeft.size() == rRight.size();
    if (output) {
        auto itLeft = rLeft.begin();
        auto itRight = rRight.begin();
        const auto itRightEnd = rRight.end();
        for (; itRight!=itRightEnd; ++itLeft, ++itRight) {
            if (!_componentComparison.equal(*itLeft, *itRight)) {
                output = false;
                break;
            }
        }
    }
    return output;
}


template <class TContainer>
requires concepts::Container<TContainer>
bool Comparison<TContainer>::less(Ref<const TContainer> rLeft,
                                  Ref<const TContainer> rRight) const noexcept
{
    auto itLeft = rLeft.begin();
    auto itRight = rRight.end();
    const auto itLeftEnd = rLeft.end();
    const auto itRightEnd = rRight.end();
    while (true) {
        if (itLeft == itLeftEnd) {
            return itRight != itRightEnd;
        } else if (itRight == itRightEnd) {
            return false;
        } else if (_componentComparison.less(*itLeft++, *itRight++)) {
            return true;
        }
    }
}


template <class T>
template <class ...TArgs>
EqualityComparison<T>::EqualityComparison(TArgs&&... rArgs) noexcept
    : _wrapped(std::forward<TArgs>(rArgs)...)
{
}


template <class T>
bool EqualityComparison<T>::operator()(Ref<const T> rLeft,
                                       Ref<const T> rRight) const noexcept
{
    return _wrapped.equal(rLeft, rRight);
}


template <class T>
template <class ...TArgs>
Ordering<T>::Ordering(TArgs&&... rArgs) noexcept
    : _wrapped(std::forward<TArgs>(rArgs)...)
{
}


template <class T>
bool Ordering<T>::operator()(Ref<const T> rLeft,
                             Ref<const T> rRight) const noexcept
{
    return _wrapped.less(rLeft, rRight);
}


} // namespace cie::utils


#endif
