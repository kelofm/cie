#ifndef CIE_UTILS_COMPARISON_HPP
#define CIE_UTILS_COMPARISON_HPP

// --- Utility Includes ---
#include "packages/compile_time/packages/concepts/inc/container_concepts.hpp"

// --- STL Includes ---
#include <concepts>


namespace cie::utils {


/// @addtogroup cieutils
template <class TInteger>
struct IntegerComparison
{
    bool equal(const TInteger left,
               const TInteger right) const noexcept;

    bool less(const TInteger left,
              const TInteger right) const noexcept;
}; // struct IntegerComparison


/// @addtogroup cieutils
template <class TFloat>
class FloatComparison
{
public:
    FloatComparison() noexcept;

    FloatComparison(const TFloat absoluteTolerance,
                    const TFloat relativeTolerance) noexcept;

    bool equal(const TFloat left,
               const TFloat right) const noexcept;

    bool less(const TFloat left,
              const TFloat right) const noexcept;

private:
    TFloat _absoluteTolerance;

    TFloat _relativeTolerance;
}; // class Comparison


/// @addtogroup cieutils
template <class T>
class Comparison {};


/// @addtogroup cieutils
template <std::floating_point TFloat>
class Comparison<TFloat> : public FloatComparison<TFloat>
{
public:
    using FloatComparison<TFloat>::FloatComparison;
};


/// @addtogroup cieutils
template <class TContainer>
requires concepts::Container<TContainer>
class Comparison<TContainer>
{
public:
    template <class ...TArgs>
    Comparison(TArgs&&... rArgs) noexcept;

    /// @brief Componentwise comparison.
    bool equal(Ref<const TContainer> rLeft,
               Ref<const TContainer> rRight) const noexcept;

    /// @brief Lexicographical comparison.
    bool less(Ref<const TContainer> rLeft,
              Ref<const TContainer> rRight) const noexcept;

private:
    Comparison<typename TContainer::value_type> _componentComparison;
}; // class Comparison<TContainer>


/// @addtogroup cieutils
template <class T>
class EqualityComparison
{
public:
    template <class ...TArgs>
    EqualityComparison(TArgs&&... rArgs) noexcept;

    bool operator()(Ref<const T> rLeft,
                    Ref<const T> rRight) const noexcept;

private:
    Comparison<T> _wrapped;
};


/// @addtogroup cieutils
template <class T>
class Ordering
{
public:
    template <class ...TArgs>
    Ordering(TArgs&&... rArgs) noexcept;

    bool operator()(Ref<const T> rLeft,
                    Ref<const T> rRight) const noexcept;

private:
    Comparison<T> _wrapped;
};


} // namespace cie::utils

#include "packages/maths/impl/Comparison_impl.hpp"

#endif
