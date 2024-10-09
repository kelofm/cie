#ifndef CIE_HIERARCHICAL_INTEGER_HPP
#define CIE_HIERARCHICAL_INTEGER_HPP

// --- Utility Includes ---
#include "packages/compile_time/packages/type_select/inc/integer.hpp"
#include "packages/compile_time/packages/concepts/inc/basic_concepts.hpp"
#include "packages/types/inc/types.hpp"

// --- STL Includes ---
#include <limits>


namespace cie {


/// @todo incomplete, tests fail - fixme
/// @ingroup cieutils
template <unsigned short Dim, concepts::UnsignedInteger TValue = Size>
class HierarchicalInteger
{
public:
    static constexpr const Size Dimension = Dim;

    using Value = TValue;

    /// @brief Shortest unsigned integer capable of representing the highest cell index derived from the highest value of @a TValue.
    using Cell = typename ct::Unsigned::ShortestWithRange<0,(std::numeric_limits<TValue>::max() >> Dim)>;

    /// @brief Shortest unsigned integer capable of representing the highest level derived from the highest value of @a TValue.
    using Level = typename ct::Unsigned::ShortestWithRange<0,(bitSize<Cell> - Dim)>;

    /// @brief Shortest unsigned integer capable of representing 2^Dimension
    using Local = typename ct::Unsigned::ShortestWithRange<0,(1 << Dim)>;

    /// @brief Maximum number of fully representable node levels.
    inline static constexpr const Level MaxLevels = bitSize<Cell> - Dimension;

public:
    /// @brief Default constructor.
    /// @note Initializes the wrapped value to 0.
    constexpr HierarchicalInteger();

    /// @brief Construct a hierarchical index from an integer.
    /// @note Triggers a decay.
    HierarchicalInteger(TValue);

    constexpr HierarchicalInteger(RightRef<HierarchicalInteger> r_rhs) = default;

    constexpr HierarchicalInteger(Ref<const HierarchicalInteger> r_rhs) = default;

    constexpr HierarchicalInteger& operator=(RightRef<HierarchicalInteger> r_rhs) = default;

    HierarchicalInteger& operator=(Ref<HierarchicalInteger> r_rhs) = default;

    /// @brief Assign a new value to the wrapped variable.
    /// @note Triggers a decay.
    HierarchicalInteger& operator=(TValue value);

    /// @brief Access the wrapped value.
    constexpr operator TValue () const noexcept;

    /// @brief Access the wrapped value.
    constexpr TValue value() const noexcept;

    /// @brief Index of the cell this node belongs to.
    /// @return Stored value / 2^Dimension
    constexpr Cell getCell() const noexcept;

    /// @brief Index of the cell @a node belongs to.
    /// @return @a node / 2^Dimension
    static constexpr Cell getCell(TValue node) noexcept;

    /// @brief Index of the cell that contains the subcell this node belongs to.
    /// @note Returns itself if invoked on the root.
    constexpr Cell getParentCell() const noexcept;

    /// @brief Index of the cell that contains the subcell of @a node.
    /// @note Returns @a node if invoked on the root.
    static constexpr Cell getParentCell(Cell cell) noexcept;

    /// @brief Find which level this node is located in.
    Level getLevel() const noexcept;

    /** @brief Find which level @a node is located in.
     *  @details Level boundaries:
     *           - Number of nodes in a cell: 2^D.
     *           - Number of cells in one level: 2^D^level.
     *           - Number of nodes in one level: 2^D^(level + 1)
     *           - Total number of nodes in a level and all levels below:
     *             000001000..01000..01000..01000..0
     *                  |     ||     ||     ||     |
     *                  |--D--||--D--||--D--||--D--|
     *                  |                          |
     *                  |-------- level * D -------|
     */
    static Level getLevel(TValue node) noexcept;

    /// @brief Get the index of this node within its containing cell.
    Local constexpr getLocalNode() const noexcept;

    /** @brief Get the index of @a node within its containing cell.
     *  @details Example layout in 2D:
     *           2 --- 3
     *           |     |
     *           0 --- 1
     */
    static constexpr Local getLocalNode(TValue node) noexcept;

    /// @brief Get the index of this node's cell within its parent.
    constexpr Local getLocalCell() const noexcept;

    /** @brief Get the index of @a node 's cell within its parent.
     *  @details Example layout in 2D:
     *           + --- + --- +
     *           |  2  |  3  |
     *           + --- + --- +
     *           |  0  |  1  |
     *           + --- + --- +
     */
    static constexpr Local getLocalCell(Cell cell) noexcept;

    /** @brief Convert the input integer to its level-wise local representation.
     *  @param node node to get the level-wise representation of.
     *  @param level level of @p node.
     *  @note The level-wise local representation is ambiguous for fixed-size integers,
     *        which is why the precomputed level is required as an input.
     */
    static TValue getLocal(TValue node, Level level) noexcept;

    static TValue fromLocal(TValue local, Level level) noexcept;

    static constexpr Cell numberOfCellsInLevel(Level level) noexcept;

    static constexpr TValue numberOfNodesInLevel(Level level) noexcept;

    static constexpr Cell numberOfCellsUpToLevel(Level level) noexcept;

    static constexpr TValue numberOfNodesUpToLevel(Level level) noexcept;

private:
    /** @brief Reduce the integer to a coincident node within the context of its parent cell.
     *  @details Example in 2D:
     *           The figure below shows a cell with one sublevel. Subcells are indexed locally,
     *           and nodes are in turn indexed locally within them. Each node is uniquely identified
     *           by its containing cell's local index and its own index within. All indices are
     *           represented in binary.
     *           10 ------------ + ------------ 11                      Decay Table
     *           | 10         11 | 10         11 |          + ----------- + --- + ----------- +
     *           |               |               |          |    From     | ==> |     To      |
     *           |      10       |      11       |          | ----------- + --- + ----------- +
     *           |               |               |          | Cell | Node |     | Cell | Node |
     *           | 00         01 | 00         01 |          | ---- + ---- + --- + ---- + ---- +
     *           + ------------- + ------------- +          |  01     00  |     |  00     01  |
     *           | 10         11 | 10         11 |          |  10     00  |     |  00     10  |
     *           |               |               |          |  10     01  |     |  00     11  |
     *           |      00       |      01       |          |  11     00  |     |  00     11  |
     *           |               |               |          |  11     01  |     |  01     11  |
     *           | 00         01 | 00         01 |          |  11     10  |     |  10     11  |
     *           00 ------------ + ------------ 01          |  00     00  |     | root    00  |
     *                                                      |  01     01  |     | root    01  |
     *                                                      |  10     10  |     | root    10  |
     *                                                      |  11     11  |     | root    11  |
     *                                                      + ---- + ---- + --- + ---- + ---- +
     */
    bool decayOnce();

    /// @brief Reduce the integer to coincident nodes until the one with the smallest index is found.
    void decay();

private:
    TValue _value;

    /// @brief Mask for simplified mod operations
    /// @note (@a value % 2^Dimension) is equivalent to (@a value & &a modMask)
    template <concepts::UnsignedInteger TI>
    inline static constexpr const TI modMask = ~(~0u << Dim);
}; // class HierarchicalInteger


} // namespace cie

#include "packages/trees/impl/HierarhicalInteger_impl.hpp"

#endif
