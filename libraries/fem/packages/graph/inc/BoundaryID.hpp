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
     *                              (1, true)
     *
     *                                  y
     *                                  ^
     *                                  |
     *                           + ---- 1 ---- +
     *                           |      |      |
     *                           |      |      |
     *           (0, false)  -- -1 ---- + ---- 1 -- >x   (0, true)
     *                           |      |      |
     *                           |      |      |
     *                           + --- -1 ---- +
     *                                  |
     *                                  |
     *
     *                              (1, false)
     *           @endcode
     *           In this example, the order of boundaries is:
     *           @code
     *           0: (0, false)
     *           1: (0,  true)
     *           2: (1, false)
     *           3: (1,  true)
     *           @endcode
     *           In this example, the order of boundaries is:
     *           @code
     *           0: (0, false)
     *           1: (0,  true)
     *           2: (1, false)
     *           3: (1,  true)
     *           @endcode
     */
    constexpr BoundaryID(unsigned dimension, bool direction);

    /** @brief Construct from dimension index and direction.
     *  @param name full name of the boundary consisting of the:
     *              - direction ("-" or "+")
     *              - dimension ("x", "y", or "z")
     *              The current implementation only allows the following names:
     *              - "-x"
     *              - "+x"
     *              - "-y"
     *              - "+y"
     *              - "-z"
     *              - "+z"
     *  @details Example in 2D:
     *           @code
     *                          "+y"
     *
     *                            y
     *                            ^
     *                            |
     *                     + ---- 1 ---- +
     *                     |      |      |
     *                     |      |      |
     *           "-x"  -- -1 ---- + ---- 1 -- >x   "+x"
     *                     |      |      |
     *                     |      |      |
     *                     + --- -1 ---- +
     *                            |
     *                            |
     *
     *                          "-y"
     *           @endcode
     *           In this example, the order of boundaries is:
     *           @code
     *           0: "-x"
     *           1: "+x"
     *           2: "-y"
     *           3: "+y"
     *           @endcode
     */
    constexpr BoundaryID(const char name[3]);

    /// @brief Increment to the next boundary.
    /// @details The next boundary is either the positive side of the
    ///          current dimension or the negative side of the next.
    constexpr BoundaryID& operator++() noexcept;

    /// @brief Increment to the next boundary.
    /// @details The next boundary is either the positive side of the
    ///          current dimension or the negative side of the next.
    constexpr BoundaryID operator++(int) noexcept;

    /// @brief Flip the boundary to represent the opposite one in the same direction.
    constexpr void flip() noexcept;

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

    /// @brief Make @ref BoundaryID sortable.
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
