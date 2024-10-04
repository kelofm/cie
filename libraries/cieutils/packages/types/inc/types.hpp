#ifndef CIE_CIEUTILS_TYPES_HPP
#define CIE_CIEUTILS_TYPES_HPP

// --- Utility Includes ---
#include "packages/macros/inc/debug.hpp"

// --- STL Includes ---
#include <cstddef>
#include <string>
#include <type_traits>
#include <climits>


namespace cie {


/// @ingroup cieutils
template <class T>
class NoInitWrapper
{
public:
    using Value = T;

public:
    NoInitWrapper() noexcept {};

    NoInitWrapper(T value) noexcept
        : _value(value)
    {}

    NoInitWrapper(NoInitWrapper&& r_rhs) noexcept = default;

    NoInitWrapper(const NoInitWrapper& r_rhs) noexcept = default;

    NoInitWrapper& operator=(T r_rhs) noexcept
    {_value = r_rhs; return *this;}

    NoInitWrapper& operator=(NoInitWrapper&& r_rhs) noexcept = default;

    NoInitWrapper& operator=(const NoInitWrapper& r_rhs) noexcept = default;

    operator T() const
    {return _value;}

    operator T&()
    {return _value;}

private:
    T _value;
}; // struct Wrapper


using Bool = bool;

using Size = size_t;

using Float = float;

using Double = double;

using String = std::string;

template <class T>
using Ptr = T*;

template <class T>
using Ref = T&;

template <class T>
using RightRef = T&&;

template <class T>
inline constexpr Size bitSize = CHAR_BIT * sizeof(T);

template <class T>
Ptr<const T> makePtrTo(Ref<const T> r_t);

template <class T>
Ptr<T> makePtrTo(Ref<T> r_t);

template <class T>
struct VoidSafe
{
    using RightRef = T&&;

    using Ref = T&;

    using Ptr = T*;
};

template <class TVoid>
requires std::is_same_v<std::remove_const_t<TVoid>,void>
struct VoidSafe<TVoid>
{
    using RightRef = void;

    using Ref = void;

    using Ptr = void;
};

} // namespace cie

#include "packages/types/impl/types_impl.hpp"
#include "packages/macros/inc/logging.hpp"

#endif
