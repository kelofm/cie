#ifndef CIE_HIERARCHICAL_INTEGER_IMPL_HPP
#define CIE_HIERARCHICAL_INTEGER_IMPL_HPP

// --- Utility Includes ---
#include "packages/trees/inc/HierarchicalInteger.hpp"
#include "packages/types/inc/types.hpp"
#include "packages/maths/inc/bit.hpp"

// --- STL Includes ---
#include <algorithm>


namespace cie {


template <unsigned short D, concepts::UnsignedInteger TValue>
inline constexpr HierarchicalInteger<D,TValue>::HierarchicalInteger()
    : _value(0)
{
}


template <unsigned short D, concepts::UnsignedInteger TValue>
inline HierarchicalInteger<D,TValue>::HierarchicalInteger(TValue value)
    : _value(value)
{
    this->decay();
}


template <unsigned short D, concepts::UnsignedInteger TValue>
inline HierarchicalInteger<D,TValue>&
HierarchicalInteger<D,TValue>::operator=(TValue value)
{
    _value = value;
    this->decay();
    return *this;
}


template <unsigned short D, concepts::UnsignedInteger TValue>
inline constexpr
HierarchicalInteger<D,TValue>::operator TValue() const noexcept
{
    return _value;
}


template <unsigned short D, concepts::UnsignedInteger TValue>
inline constexpr TValue
HierarchicalInteger<D,TValue>::value() const noexcept
{
    return _value;
}


template <unsigned short D, concepts::UnsignedInteger TValue>
inline constexpr typename HierarchicalInteger<D,TValue>::Cell
HierarchicalInteger<D,TValue>::getCell() const noexcept
{
    return HierarchicalInteger::getCell(_value);
}


template <unsigned short D, concepts::UnsignedInteger TValue>
inline constexpr typename HierarchicalInteger<D,TValue>::Cell
HierarchicalInteger<D,TValue>::getCell(TValue node) noexcept
{
    return node >> D;
}


template <unsigned short D, concepts::UnsignedInteger TValue>
inline constexpr typename HierarchicalInteger<D,TValue>::Cell
HierarchicalInteger<D,TValue>::getParentCell() const noexcept
{
    return HierarchicalInteger::getParentCell(this->getCell());
}


template <unsigned short D, concepts::UnsignedInteger TValue>
inline constexpr typename HierarchicalInteger<D,TValue>::Cell
HierarchicalInteger<D,TValue>::getParentCell(Cell cell) noexcept
{
    // Handle integer underflow
    return ((cell != 0u) * (cell - 1u)) >> D;
}


template <unsigned short D, concepts::UnsignedInteger TValue>
inline typename HierarchicalInteger<D,TValue>::Level
HierarchicalInteger<D,TValue>::getLevel() const noexcept
{
    return HierarchicalInteger::getLevel(_value);
}


template <unsigned short D, concepts::UnsignedInteger TValue>
inline typename HierarchicalInteger<D,TValue>::Level
HierarchicalInteger<D,TValue>::getLevel(TValue node) noexcept
{
    static_assert(D <= bitSize<TValue>);

    const TValue cellIndex = node >> D;
    TValue cellsInLevel = 0;
    Level level = 0;

    while ((cellsInLevel |= (cellsInLevel << D) | 1u) <= cellIndex && level < HierarchicalInteger::MaxLevels) ++level;
    return level;
}


template <unsigned short D, concepts::UnsignedInteger TValue>
inline constexpr typename HierarchicalInteger<D,TValue>::Local
HierarchicalInteger<D,TValue>::getLocalNode() const noexcept
{
    return HierarchicalInteger::getLocalNode(_value);
}


template <unsigned short D, concepts::UnsignedInteger TValue>
inline constexpr typename HierarchicalInteger<D,TValue>::Local
HierarchicalInteger<D,TValue>::getLocalNode(TValue node) noexcept
{
    // Equivalent to: node % 2^D
    return node & HierarchicalInteger::modMask<TValue>;
}


template <unsigned short D, concepts::UnsignedInteger TValue>
inline constexpr typename HierarchicalInteger<D,TValue>::Local
HierarchicalInteger<D,TValue>::getLocalCell() const noexcept
{
    return HierarchicalInteger::getLocalCell(this->getCell());
}


template <unsigned short D, concepts::UnsignedInteger TValue>
inline constexpr typename HierarchicalInteger<D,TValue>::Local
HierarchicalInteger<D,TValue>::getLocalCell(Cell cell) noexcept
{
    return ((cell | (cell == 0)) - 1u) & HierarchicalInteger::modMask<Cell>; // avoid underflow
}


template <unsigned short D, concepts::UnsignedInteger TValue>
inline TValue
HierarchicalInteger<D,TValue>::getLocal(TValue node, Level level) noexcept
{
    TValue local = HierarchicalInteger::getLocalNode(node);
    Cell cell = HierarchicalInteger::getCell(node);
    for (Level i_level=1; i_level<level; ++i_level, cell=HierarchicalInteger::getParentCell(cell))
        local |= (HierarchicalInteger::getLocalCell(cell)) << (D * i_level);
    return local;
}


template <unsigned short D, concepts::UnsignedInteger TValue>
inline TValue
HierarchicalInteger<D,TValue>::fromLocal(TValue local, Level level) noexcept
{
    // Level cannot be 0 (no node is on cell level 0)
    //assert(level);
    --level;
    TValue offset = level * D;
    TValue node = (local >> offset) + (local && level);
    level -= bool(level);
    offset -= D * bool(offset);

    TValue mask = ~(~0u << D) << offset;

    if (offset)
        for(; mask; mask >>= D, offset-=D)
        {
            node <<= D;
            node |= (mask & local) >> offset;
        }

    return node;
}


template <unsigned short D, concepts::UnsignedInteger TValue>
inline constexpr typename HierarchicalInteger<D,TValue>::Cell
HierarchicalInteger<D,TValue>::numberOfCellsInLevel(Level level) noexcept
{
    //assert(level + D < bitSize<Cell>);
    return (Cell(1) << D) << level;
}


template <unsigned short D, concepts::UnsignedInteger TValue>
inline constexpr TValue
HierarchicalInteger<D,TValue>::numberOfNodesInLevel(Level level) noexcept
{
    return TValue(HierarchicalInteger<D,TValue>::numberOfCellsInLevel(level)) << D;
}


template <unsigned short D, concepts::UnsignedInteger TValue>
inline constexpr typename HierarchicalInteger<D,TValue>::Cell
HierarchicalInteger<D,TValue>::numberOfCellsUpToLevel(Level level) noexcept
{
    Cell cells = 1;
    for (Level l=0; l<level; ++l) cells |= (cells << D);
    return cells;
}


template <unsigned short D, concepts::UnsignedInteger TValue>
inline constexpr TValue
HierarchicalInteger<D,TValue>::numberOfNodesUpToLevel(Level level) noexcept
{
    return TValue(HierarchicalInteger<D,TValue>::numberOfCellsUpToLevel(level)) << D;
}


template <unsigned short D, concepts::UnsignedInteger TValue>
inline bool
HierarchicalInteger<D,TValue>::decayOnce()
{
    const TValue old = _value;
    const Cell cell = this->getCell();
    const Local localCell = HierarchicalInteger::getLocalCell(cell);
    const Local localNode = this->getLocalNode();

    {
        // Local cell and node indices line up
        // => node can decay to the parent cell's node
        //    at the same position.
        //_value = (TValue(this->getParentCell()) << D) + localNode;

        // Disassemble node
        Level level = this->getLevel() + 1;
        const TValue local = HierarchicalInteger::getLocal(_value, level);
        TValue newLocal = local;

        // Find the first occasion where the following condition is true:
        // - the parent's bit is set but the child's is not.
        for (unsigned short d=0; d<D && newLocal==local; ++d)
        {
            TValue parentMask = 1 << ((level + 1) * D - d - 1);
            TValue childMask = parentMask >> D;
            do {
                Local pm = parentMask & newLocal;
                Local cm = (~(childMask & newLocal)) << D;
                newLocal ^= pm & cm;
            } while ((parentMask>>=D) | (childMask>>=D));
        }

        // Reassemble node
        _value = this->fromLocal(newLocal, level);

        return _value != old;
    }

    if (localNode < localCell)
    {
        // Decay the node within its cell's parent.

        // Equivalent algorithm:
        // Loop through the (binary) digits of the local cell and node.
        // If the cell's digit is greater than that of the node, swap them.
        Local mask = localCell ^ localNode & localCell;

        // Recompute the node
        const Cell newCell = this->getCell() + (localCell ^ mask) - localCell;
        _value = (newCell << D) + (localNode ^ mask);
    }

    return _value != old;
}


template <unsigned short D, concepts::UnsignedInteger TValue>
inline void
HierarchicalInteger<D,TValue>::decay()
{
    while (this->decayOnce()) {}
}


} // namespace cie


#endif
