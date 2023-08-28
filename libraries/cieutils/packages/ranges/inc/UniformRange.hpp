#ifndef CIE_UTILS_RANGES_UNIFORM_RANGE_HPP
#define CIE_UTILS_RANGES_UNIFORM_RANGE_HPP

// --- Utility Includes ---
#include "packages/compile_time/packages/concepts/inc/basic_concepts.hpp"
#include "packages/macros/inc/checks.hpp"
#include "packages/macros/inc/typedefs.hpp"


namespace cie::utils {


template <concepts::Numeric T>
class UniformRange
{
private:
    class Iterator
    {
    public:
        using value_type = T;

        using difference_type = std::ptrdiff_t;

        using reference = T&;

        using pointer = T*;

        using iterator_category = std::random_access_iterator_tag;

    public:
        Iterator()
            : Iterator(nullptr)
        {CIE_THROW(NullPtrException, "")}

        Iterator(const UniformRange* p_range)
            : Iterator(p_range->_begin, p_range)
        {}

        Iterator(Size state, const UniformRange* p_range)
            : _state(state),
              _p_range(p_range)
        {}

        Iterator(Iterator&& r_rhs) = default;

        Iterator(const Iterator& r_rhs) = default;

        Iterator& operator=(Iterator&& r_rhs) = default;

        Iterator& operator=(const Iterator& r_rhs) = default;

        value_type operator*()
        {return _p_range->_begin + T(_state) * _p_range->_step;}

        value_type operator*() const
        {return _p_range->_begin + T(_state) * _p_range->_step;}

        Iterator& operator++()
        {++_state; return *this;}

        Iterator operator++(int)
        {Iterator copy = *this; ++(*this); return copy;}

        Iterator& operator--()
        {--_state; return *this;}

        Iterator operator--(int)
        {Iterator copy = *this; --(*this); return copy;}

        Iterator& operator+=(difference_type difference)
        {_state += difference; return *this;}

        Iterator& operator-=(difference_type difference)
        {_state -= difference; return *this;}

        friend Iterator operator+(const Iterator& r_lhs, difference_type rhs)
        {return Iterator(r_lhs._state + rhs, r_lhs._p_range);}

        friend Iterator operator+(difference_type lhs, const Iterator& r_rhs)
        {return r_rhs + lhs;}

        friend Iterator operator-(const Iterator& r_lhs, difference_type rhs)
        {return Iterator(r_lhs._state - rhs, r_lhs._p_range);}

        friend Iterator operator-(const Iterator& r_lhs, const Iterator& r_rhs)
        {
            CIE_OUT_OF_RANGE_CHECK(r_lhs._p_range == r_rhs._p_range)
            return difference_type(r_lhs._state) - difference_type(r_rhs._state);
        }

        friend bool operator==(const Iterator& r_lhs, const Iterator& r_rhs)
        {return r_lhs._state == r_rhs._state && r_lhs._p_range == r_rhs._p_range;}

        friend bool operator!=(const Iterator& r_lhs, const Iterator& r_rhs)
        {return !(r_lhs == r_rhs);}

        friend bool operator<(const Iterator& r_lhs, const Iterator& r_rhs)
        {
            CIE_OUT_OF_RANGE_CHECK(r_lhs._p_range == r_rhs._p_range)
            return r_lhs._state < r_rhs._state;
        }

    private:
        Size _state;

        const UniformRange* _p_range;
    }; // class Iterator

public:
    using value_type = T;

    using size_type = Size;

    using iterator = Iterator;

    using const_iterator = Iterator;

    CIE_DEFINE_CLASS_POINTERS(UniformRange)

public:
    UniformRange()
        : UniformRange(T(0), T(0))
    {}

    UniformRange(T end)
        : UniformRange(T(0), end)
    {}

    UniformRange(T begin, T end)
        : UniformRange(begin, end, T(1))
    {}

    UniformRange(T begin, T end, T step)
        : _begin(begin),
          _end(end),
          _step(step)
    {
        CIE_CHECK(step != 0, "Step size cannot be 0")
        CIE_CHECK(begin < end == 0 < step, "Step sign must match the range's direction")
    }

    UniformRange(UniformRange&& r_rhs) = default;

    UniformRange(const UniformRange& r_rhs) = default;

    UniformRange& operator=(UniformRange&& r_rhs) = default;

    UniformRange& operator=(const UniformRange& r_rhs) = default;

    size_type size() const
    {return end() - begin();}

    iterator begin()
    {return iterator(0, this);}

    const iterator begin() const
    {return const_iterator(0, this);}

    iterator end()
    {return iterator(Size((_begin < _end ? _end - _begin : _begin - _end) / _step + 1), this);}

    const_iterator end() const
    {return const_iterator(Size((_begin < _end ? _end - _begin : _begin - _end) / _step + 1), this);}

private:
    T _begin;

    T _end;

    T _step;
}; // class UniformRange


} // namespace cie::utils


#endif