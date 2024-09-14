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

    StrongTypeDef(T&& r_rhs) noexcept requires concepts::MoveConstructible<T>
        : T(std::move(r_rhs))
    {}

    /// Allow move constructor if possible
    StrongTypeDef(StrongTypeDef&& r_rhs) noexcept requires concepts::MoveConstructible<T>
        : T(std::move(r_rhs))
    {}

    StrongTypeDef(const T& r_rhs) requires concepts::CopyConstructible<T>
        : T(r_rhs)
    {}

    StrongTypeDef(T& r_rhs) requires concepts::CopyConstructible<T>
        : T(r_rhs)
    {}

    /// Allow copy constructor if possible
    StrongTypeDef(const StrongTypeDef& r_rhs) requires concepts::CopyConstructible<T>
        : T(r_rhs)
    {}

    StrongTypeDef& operator=(T&& r_rhs) noexcept requires concepts::MoveAssignable<T>
    {static_cast<T&>(*this) = std::move(r_rhs); return *this;}

    /// Allow move assignment operator if possible
    StrongTypeDef& operator=(StrongTypeDef&& r_rhs) noexcept requires concepts::MoveAssignable<T>
    {static_cast<T&>(*this) = std::move(static_cast<T&&>(r_rhs)); return *this;}

    StrongTypeDef& operator=(const T& r_rhs) requires concepts::CopyAssignable<T>
    {static_cast<T&>(*this) = r_rhs; return *this;}

    StrongTypeDef& operator=(T& r_rhs) requires concepts::CopyAssignable<T>
    {static_cast<T&>(*this) = r_rhs; return *this;}

    /// Allow copy assignment operator if possible
    StrongTypeDef& operator=(const StrongTypeDef& r_rhs) requires concepts::CopyAssignable<T>
    {static_cast<T&>(*this) = static_cast<const T&>(r_rhs); return *this;}

    StrongTypeDef& operator=(StrongTypeDef& r_rhs) requires concepts::CopyAssignable<T>
    {static_cast<T&>(*this) = static_cast<const T&>(r_rhs); return *this;}

    /// Delete every assignment operator that's not in the base class, move or copy assignment operator
    template <class TT>
    StrongTypeDef& operator=(TT&& r_rhs) = delete;

    operator T&()
    {return static_cast<T&>(*this);}

    operator const T&() const
    {return static_cast<const T&>(*this);}
};


template <concepts::Integral T, class Tag>
class StrongTypeDef<T,Tag>
{
public:
    using Value = T;

    StrongTypeDef() noexcept
    {}

    StrongTypeDef(T wrapped) noexcept
        : _wrapped(wrapped)
    {}

    template <class TT>
    requires std::convertible_to<TT,T>
    StrongTypeDef(TT wrapped)
        : _wrapped(wrapped)
    {}

    StrongTypeDef& operator=(T rhs) noexcept
    {_wrapped = rhs; return *this;}

    operator const T&() const noexcept
    {return _wrapped;}

    operator T&() noexcept
    {return _wrapped;}

    template <class TT>
    requires std::convertible_to<T,TT>
    explicit operator TT() const noexcept
    {return _wrapped;}

    friend void swap(Ref<StrongTypeDef> rLeft,
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
