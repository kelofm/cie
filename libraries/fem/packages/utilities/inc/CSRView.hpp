#pragma once

// --- STL Includes ---
#include <span>
#include <type_traits>


namespace cie::fem {


template <
    class TValue,
    class TIndex = std::conditional_t<
        std::is_const_v<TValue>,
        const int,
        int>>
class CSRView {
public:
    using Value = std::remove_const_t<TValue>;

    using Index = std::remove_const_t<TIndex>;

    constexpr CSRView() noexcept = default;

    constexpr CSRView(
        std::size_t columnCount,
        std::span<TIndex> rowExtents,
        std::span<TIndex> columnIndices,
        std::span<TValue> entries) noexcept
        : _columnCount(columnCount),
          _rowExtents(rowExtents),
          _columnIndices(columnIndices),
          _entries(entries)
    {}

    std::size_t rowCount() const noexcept {
        return _rowExtents.size() - 1;
    }

    constexpr std::size_t columnCount() const noexcept {
        return _columnCount;
    }

    std::span<const Index> rowExtents() const noexcept {
        return {_rowExtents.data(), _rowExtents.size()};
    }

    std::span<Index> rowExtents() noexcept
    requires (!std::is_const_v<TIndex>) {
        return {_rowExtents.data(), _rowExtents.size()};
    }

    std::span<const Index> columnIndices() const noexcept {
        return {_columnIndices.data(), _columnIndices.size()};
    }

    std::span<Index> columnIndices() noexcept
    requires (!std::is_const_v<TIndex>) {
        return {_columnIndices.data(), _columnIndices.size()};
    }

    std::span<const Value> entries() const noexcept {
        return {_entries.data(), _entries.size()};
    }

    std::span<Value> entries() noexcept
    requires (!std::is_const_v<TValue>) {
        return {_entries.data(), _entries.size()};
    }

private:
    std::size_t _columnCount;

    std::span<TIndex> _rowExtents, _columnIndices;

    std::span<TValue> _entries;
}; // struct CSRView


} // namespace cie::fem
