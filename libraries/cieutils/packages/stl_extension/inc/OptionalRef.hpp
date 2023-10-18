#ifndef CIE_UTILS_OPTIONAL_REF_HPP
#define CIE_UTILS_OPTIONAL_REF_HPP

// --- Utility Includes ---
#include "packages/macros/inc/exceptions.hpp"

// --- STL Includes ---
#include <optional>
#include <stdexcept>
#include <type_traits> // alignment_of
#include <bitset> // bitset
#include <cstdint> // uintptr_t


namespace cie {


/// @brief An optional, modeling a pointer.
/// @details @a OptionalRef stores a pointer to the instance it was assigned
///          and stores the value indicator in the least significant bit of
///          the pointer. As a result, the memory alignment of @a OptionalRef
///          is identical to that of the stored pointer.
/// @addtogroup cieutils
template <class TElement>
class OptionalRef
{
private:
    using PointerInteger = std::uintptr_t;

public:
    using value_type = TElement&;

    /// @brief Construct an invalid (uninitialized) optional.
    OptionalRef() noexcept
    {
        // First of all, make sure that the pointer doesn't use the least significant bit
        static_assert(1 < std::alignment_of<TElement>(),
                      "OptionalReference cannot be used with types of size 1");

        // Assume that the size of PointerInteger is equal to the system's pointer size
        static_assert(sizeof(PointerInteger) == sizeof(TElement*));

        // Verify that the alignment is identical to a pointer
        static_assert(std::alignment_of<OptionalRef>() == std::alignment_of<TElement*>());

        mpElement.pointer = 0;
    }

    /// @brief Construct an optional reference to the provided instance.
    OptionalRef(TElement& rElement) noexcept
    {
        PointerInteger p = reinterpret_cast<PointerInteger>(&rElement);
        mpElement.integer = p | static_cast<PointerInteger>(1);
    }

    OptionalRef(OptionalRef&& rRhs) noexcept = default;

    OptionalRef(const OptionalRef& rRhs) noexcept = default;

    OptionalRef& operator=(OptionalRef&& rRhs) noexcept = default;

    OptionalRef& operator=(const OptionalRef& rRhs) noexcept = default;

    /// @brief Assign the provided instance to the optional.
    OptionalRef& operator=(TElement& rElement) noexcept
    {
        PointerInteger p = reinterpret_cast<PointerInteger>(&rElement);
        mpElement.integer = p | static_cast<PointerInteger>(1);
    }

    /// @brief Check whether the optional is valid.
    bool has_value() const noexcept
    {
        return mpElement.integer & static_cast<PointerInteger>(1);
    }

    /// @brief Check whether the optional is valid.
    explicit operator bool () const noexcept
    {
        return this->has_value();
    }

    /// @brief Access the stored reference.
    /// @throws if the optional is uninitialized.
    TElement& value() const
    {
        if (!this->has_value()){
            CIE_THROW(Exception, "Bad optional access")
        } else {
            PointerInteger p = mpElement.integer & ~static_cast<PointerInteger>(1);
            return *reinterpret_cast<TElement*>(p);
        }
    }

    /// @brief Access the stored reference.
    /// @throws if the optional is uninitialized.
    TElement& operator*() const noexcept
    {
        return this->value();
    }

    /// @brief Clear the optional.
    void reset() noexcept
    {
        mpElement.integer = 0;
    }

    /// @brief Cast to constant reference
    operator OptionalRef<const TElement> () noexcept
    {
        if (this->has_value()) {
            return OptionalRef<const TElement>(this->value());
        } else {
            return {};
        }
    }

private:
    union {
        TElement* pointer;
        PointerInteger integer;
    } mpElement;
}; // class OptionalRef


} // namespace cie

#endif
