#ifndef CIE_FEM_GRAPH_BOUNDARY_ID_HPP
#define CIE_FEM_GRAPH_BOUNDARY_ID_HPP

// --- FEM Includes ---
#include "packages/io/inc/GraphML.hpp" // GraphML::Deserializer

// --- Utility Includes ---
#include "packages/types/inc/types.hpp"

// --- STL Includes ---
#include <ostream> // std::ostream
#include <functional> // std::hash
#include <cstdint> // std::uint8_t


namespace cie::fem {


/// @brief Class for identifying hypercubes' faces (left, right, bottom, top, etc.) in local space.
class BoundaryID {
public:
    // Required by std::incrementable.
    using difference_type = std::int8_t;

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
    [[nodiscard]] constexpr BoundaryID operator-() const noexcept;

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
    std::uint8_t _id;
}; // class BoundaryID


/// @brief Get a string representation of the boundary ("+0", "-1", ...).
Ref<std::ostream> operator<<(Ref<std::ostream> rStream, BoundaryID id);


inline bool operator==(std::pair<cie::fem::BoundaryID,cie::fem::BoundaryID> left,
                       std::pair<cie::fem::BoundaryID,cie::fem::BoundaryID> right) noexcept
{
    return (left.first == right.first) && (left.second == right.second);
}


inline bool operator!=(std::pair<cie::fem::BoundaryID,cie::fem::BoundaryID> left,
                       std::pair<cie::fem::BoundaryID,cie::fem::BoundaryID> right) noexcept
{
    return left.first != right.first && left.second != right.second;
}


} // namespace cie::fem



// --- IO --- //



namespace cie::fem::io {


template <>
struct io::GraphML::Serializer<BoundaryID>
{
    void header(Ref<XMLElement> rElement);

    void operator()(Ref<XMLElement> rElement,
                    Ref<const BoundaryID> rInstance);
}; // struct GraphML::Serializer<BoundaryID>


template <>
struct io::GraphML::Deserializer<BoundaryID>
    : public io::GraphML::DeserializerBase<BoundaryID>
{
    using io::GraphML::DeserializerBase<BoundaryID>::DeserializerBase;

    static void onElementBegin(Ptr<void> pThis,
                               std::string_view elementName,
                               std::span<GraphML::AttributePair> attributes) noexcept;

    static void onText(Ptr<void> pThis,
                       std::string_view data);

    static void onElementEnd(Ptr<void> pThis,
                             std::string_view elementName) noexcept;
}; // struct GraphML::Deserializer<BoundaryID>


} // namespace cie::fem::io



// --- Hash --- //



namespace std {


/// @brief Standard @a hash specialization for @ref cie::fem::BoundaryID.
template <>
struct hash<cie::fem::BoundaryID>
{
    size_t operator()(cie::fem::BoundaryID id) const noexcept;
}; // hash<BoundaryID>


template <>
struct hash<pair<cie::fem::BoundaryID,cie::fem::BoundaryID>>
{
    auto operator()(pair<cie::fem::BoundaryID,cie::Size> item) const
    {
        const auto tmp = hash<cie::fem::BoundaryID>()(item.first);
        return tmp ^ (hash<cie::Size>()(item.second) + 0x9e3779b9 + (tmp<<6) + (tmp>>2)); // <== from boost::hash_combine
    }
}; // hash<pair<BoundaryID,BoundaryID>>


} // namespace std

#include "packages/graph/impl/BoundaryID_impl.hpp"

#endif
