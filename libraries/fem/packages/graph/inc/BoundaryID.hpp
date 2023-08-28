#ifndef CIE_FEM_GRAPH_BOUNDARY_ID_HPP
#define CIE_FEM_GRAPH_BOUNDARY_ID_HPP

// --- Utility Includes ---
#include "packages/macros/inc/checks.hpp"
#include "packages/maths/inc/bit.hpp"
#include "packages/types/inc/types.hpp"

// --- STL Includes ---
#include <ostream> // std::ostream
#include <functional> // std::hash


namespace cie::fem {


/// @brief Class for identifying hypercubes' faces (left, right, bottom, top, etc.) in local space.
class BoundaryID
{
public:
    /// @brief Construct an ID for the negative side (-1) of the first dimension (x).
    constexpr BoundaryID() noexcept;

    /** @brief Construct from dimension index and direction.
     *  @param dimension Index of the dimension, beginning with 0 (x:0, y:1, ...).
     *  @param direction Index of the direction (negative: false, positive: true).
     *  @details Example in 2D:
     *           @code
     *                        (1,true)
     *                      +----------+
     *                      |          |
     *           (0, false) |          | (0, true)
     *                      |          |
     *                      +----------+
     *                       (1, false)
     *           @endcode
     */
    constexpr BoundaryID(unsigned dimension, bool direction);

    /// @brief Increment to the next boundary.
    /// @details The next boundary is either the positive side of the
    ///          current dimension or the negative side of the next.
    constexpr BoundaryID& operator++() noexcept;

    /// @brief Increment to the next boundary.
    /// @details The next boundary is either the positive side of the
    ///          current dimension or the negative side of the next.
    constexpr BoundaryID operator++(int) noexcept;

    /// @brief Get the dimension index (0:x, 1:y, ...).
    unsigned getDimension() const noexcept;

    /// @brief Get the direction index (false: negative, true: positive)
    constexpr bool getDirection() const noexcept;

    /// @brief Get the underlying representation.
    constexpr unsigned getData() const noexcept;

    /// @brief Check whether two instances refer to the same boundary.
    constexpr friend bool operator==(BoundaryID left, BoundaryID right) noexcept;

    /// @brief Check whether two instances refer to separate boundaries.
    constexpr friend bool operator!=(BoundaryID left, BoundaryID right) noexcept;

    constexpr friend bool operator<(BoundaryID left, BoundaryID right) noexcept;

private:
    unsigned _id;
}; // class BoundaryID


/// @brief Get a string representation of the boundary ("+x", "-y", ...).
Ref<std::ostream> operator<<(Ref<std::ostream> r_stream, BoundaryID id);


} // namespace cie::fem


namespace std {


/// @brief Standard @a hash specialization for @ref BoundaryID.
template <>
struct hash<cie::fem::BoundaryID>
{
    size_t operator()(cie::fem::BoundaryID id) const noexcept;
}; // hash<BoundaryID>


} // namespace std

#include "packages/graph/impl/BoundaryID_impl.hpp"

#endif
