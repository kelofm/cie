#ifndef CIE_UTILS_STL_EXTENSION_STRONG_TYPEDEF_HPP
#define CIE_UTILS_STL_EXTENSION_STRONG_TYPEDEF_HPP

// --- Utility Includes ---
#include "packages/compile_time/packages/concepts/inc/basic_concepts.hpp"

// --- STL Includes ---
#include <concepts> // convertible_to
#include <iosfwd> // ostream
#include <functional> // hash


namespace cie::utils {


template <class T, class Tag>
class StrongTypeDef {};


/** Wrapper class for creating non-interchangable typedefs
 *
 *              Base
 *             /    \
 *            /      \
 *           /        \
 *          /          \
 *         /            \
 *    TypeDef1---X---TypeDef2
 *
 *  - Typedefd classes are constructible from and convertible to the Base class
 *  - Typedefd classes (TypeDef1 and TypeDef2) sharing the same Base are not
 *    implicitly convertible (though can be converted by explicitly casting through Base)
 *  - Typedefd classes share (almost*) all functionality of the Base class
 *
 *  @note * operator overloads (templated ones especially) need a bit of explicit nudging
 *  to get picked up by the compiler. Also, there's no efficient and general way to get list
 *  initialization working.
 *  @ingroup cieutils
 */
template <concepts::Deriveable T, class Tag>
class StrongTypeDef<T,Tag> : public T
{
public:
    using Value = T;

    /// Inherit all base constructors
    using T::T;

    StrongTypeDef(T&& rRhs) noexcept requires concepts::MoveConstructible<T>
        : T(std::move(rRhs))
    {}

    /// Allow move constructor if possible
    StrongTypeDef(StrongTypeDef&& rRhs) noexcept requires concepts::MoveConstructible<T>
        : T(std::move(rRhs))
    {}

    StrongTypeDef(const T& rRhs) requires concepts::CopyConstructible<T>
        : T(rRhs)
    {}

    StrongTypeDef(T& rRhs) requires concepts::CopyConstructible<T>
        : T(rRhs)
    {}

    /// Allow copy constructor if possible
    StrongTypeDef(const StrongTypeDef& rRhs) requires concepts::CopyConstructible<T>
        : T(rRhs)
    {}

    StrongTypeDef& operator=(T&& rRhs) noexcept requires concepts::MoveAssignable<T>
    {static_cast<T&>(*this) = std::move(rRhs); return *this;}

    /// Allow move assignment operator if possible
    StrongTypeDef& operator=(StrongTypeDef&& rRhs) noexcept requires concepts::MoveAssignable<T>
    {static_cast<T&>(*this) = std::move(static_cast<T&&>(rRhs)); return *this;}

    StrongTypeDef& operator=(const T& rRhs) requires concepts::CopyAssignable<T>
    {static_cast<T&>(*this) = rRhs; return *this;}

    StrongTypeDef& operator=(T& rRhs) requires concepts::CopyAssignable<T>
    {static_cast<T&>(*this) = rRhs; return *this;}

    /// Allow copy assignment operator if possible
    StrongTypeDef& operator=(const StrongTypeDef& rRhs) requires concepts::CopyAssignable<T>
    {static_cast<T&>(*this) = static_cast<const T&>(rRhs); return *this;}

    StrongTypeDef& operator=(StrongTypeDef& rRhs) requires concepts::CopyAssignable<T>
    {static_cast<T&>(*this) = static_cast<const T&>(rRhs); return *this;}

    /// Delete every assignment operator that's not in the base class, move or copy assignment operator
    template <class TT>
    StrongTypeDef& operator=(TT&& rRhs) = delete;

    operator T&()
    {return static_cast<T&>(*this);}

    operator const T&() const
    {return static_cast<const T&>(*this);}
};


template <class T, class Tag>
requires (std::is_integral_v<T> || std::is_floating_point_v<T>)
class StrongTypeDef<T,Tag>
{
public:
    using Value = T;

    constexpr StrongTypeDef() noexcept
    {}

    constexpr StrongTypeDef(T wrapped) noexcept
        : _wrapped(wrapped)
    {}

    template <class TT>
    requires std::convertible_to<TT,T>
    constexpr StrongTypeDef(TT wrapped)
        : _wrapped(wrapped)
    {}

    constexpr StrongTypeDef& operator=(T rhs) noexcept
    {_wrapped = rhs; return *this;}

    constexpr operator const T&() const noexcept
    {return _wrapped;}

    constexpr operator T&() noexcept
    {return _wrapped;}

    template <class TT>
    requires std::convertible_to<T,TT>
    explicit constexpr operator TT() const noexcept
    {return _wrapped;}

    friend constexpr void swap(Ref<StrongTypeDef> rLeft,
                               Ref<StrongTypeDef> rRight) noexcept
    {std::swap(rLeft._wrapped, rRight._wrapped);}

    friend std::ostream& operator<<(Ref<std::ostream> rStream,
                                    StrongTypeDef self)
    {return rStream << self._wrapped;}

private:
    T _wrapped;
};


} // namespace cie::utils


namespace std {


template <class T, class TTag>
struct hash<cie::utils::StrongTypeDef<T,TTag>>
{
    size_t operator()(cie::Ref<const cie::utils::StrongTypeDef<T,TTag>>& rInstance) const noexcept
    {
        using Value = typename cie::utils::StrongTypeDef<T,TTag>::Value;
        hash<Value> hasher;
        return hasher((const Value&) rInstance);
    }
}; // struct hash


} // namespace std


#define CIE_STRONG_TYPEDEF( BaseType, SubType )                                     \
    struct _typedef_ ## SubType ## _ ## BaseType ## _tag {};                        \
    typedef cie::utils::StrongTypeDef<                                              \
                            BaseType ,                                              \
                           _typedef_  ## SubType ## _ ## BaseType ## _tag> SubType


#include "packages/stl_extension/inc/StrongTypeDef_overloads.hpp"

#endif
