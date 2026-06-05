#pragma once

// --- STL Includes ---
#include <span>
#include <type_traits>


namespace cie::linalg {


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
        Index columnCount,
        std::span<TIndex> rowExtents,
        std::span<TIndex> columnIndices,
        std::span<TValue> entries) noexcept
        : _columnCount(columnCount),
          _rowExtents(rowExtents),
          _columnIndices(columnIndices),
          _entries(entries)
    {}

    TIndex rowCount() const noexcept {
        return _rowExtents.size() - 1;
    }

    constexpr TIndex columnCount() const noexcept {
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

    operator CSRView<const Value,const Index> () const noexcept
    requires (!std::is_const_v<TValue> || !std::is_const_v<TIndex>) {
        return CSRView<const Value,const Index>(
            _columnCount,
            _rowExtents,
            _columnIndices,
            _entries);
    }

    operator CSRView<Value,const Index> () noexcept
    requires (!std::is_const_v<TValue>) {
        return CSRView<Value,const Index>(
            _columnCount,
            _rowExtents,
            _columnIndices,
            _entries);
    }

private:
    TIndex _columnCount;

    std::span<TIndex> _rowExtents, _columnIndices;

    std::span<TValue> _entries;
}; // struct CSRView


template <
    class TValue,
    class TIndex = std::conditional_t<
        std::is_const_v<TValue>,
        const int,
        int>>
class CSCView : private CSRView<TValue,TIndex> {
public:
    using Value = std::remove_const_t<TValue>;

    using Index = std::remove_const_t<TIndex>;

    constexpr CSCView() noexcept = default;

    constexpr CSCView(
        TIndex rowCount,
        std::span<TIndex> columnExtents,
        std::span<TIndex> rowIndices,
        std::span<TValue> entries) noexcept
            : CSRView<TValue,TIndex>(
                rowCount,
                columnExtents,
                rowIndices,
                entries)
    {}

    TIndex rowCount() const noexcept {
        return CSRView<TValue,TIndex>::columnCount();
    }

    constexpr TIndex columnCount() const noexcept {
        return CSRView<TValue,TIndex>::rowCount();
    }

    std::span<const Index> columnExtents() const noexcept {
        return CSRView<TValue,TIndex>::rowExtents();
    }

    std::span<Index> columnExtents() noexcept
    requires (!std::is_const_v<TIndex>) {
        return CSRView<TValue,TIndex>::rowExtents();
    }

    std::span<const Index> rowIndices() const noexcept {
        return CSRView<TValue,TIndex>::columnIndices();
    }

    std::span<Index> rowIndices() const noexcept
    requires (!std::is_const_v<TIndex>) {
        return CSRView<TValue,TIndex>::columnIndices();
    }

    std::span<const Value> entries() const noexcept {
        return CSRView<TValue,TIndex>::entries();
    }

    std::span<Value> entries() noexcept
    requires (!std::is_const_v<TValue>) {
        return CSRView<TValue,TIndex>::entries();
    }

    operator CSCView<const Value,const Index>() const noexcept
    requires (!std::is_const_v<TValue> || !std::is_const_v<TIndex>) {
        return CSCView<const Value,const Index>(
            this->rowCount(),
            this->columnExtents(),
            this->rowIndices(),
            this->entries());
    }

    operator CSCView<Value,const Index>() const noexcept
    requires (!std::is_const_v<TValue>) {
        return CSCView<Value,const Index>(
            this->rowCount(),
            this->columnExtents(),
            this->rowIndices(),
            this->entries());
    }
}; // class CSCView


} // namespace cie::linalg
