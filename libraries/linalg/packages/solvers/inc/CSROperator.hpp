#pragma once

// --- Linalg Includes ---
#include "packages/solvers/inc/LinearOperator.hpp"
#include "packages/solvers/inc/DefaultSpace.hpp"

// --- Utility Includes ---
#include "packages/stl_extension/inc/OptionalRef.hpp"
#include "packages/concurrency/inc/ThreadPoolBase.hpp"

// --- STL Includes ---
#include <span>


namespace cie::linalg {


/// @brief Linear operator representing a scaled matrix-vector product in CSR format.
template <class TIndex, class TValue, class TMatrixValue = TValue>
class CSROperator
    : public LinearOperator<DefaultSpace<TValue,tags::SMP>> {
private:
    using Space = DefaultSpace<TValue,tags::SMP>;

public:
    constexpr CSROperator() noexcept = default;

    CSROperator(
        TIndex columnCount,
        std::span<const TIndex> rowExtents,
        std::span<const TIndex> columnIndices,
        std::span<const TMatrixValue> entries,
        OptionalRef<mp::ThreadPoolBase> rMaybeThreads = {});

    /// @copydoc LinearOperator::product
    void product(
        typename Space::ConstVectorView in,
        typename Space::Value scale,
        typename Space::VectorView out) override;

protected:
    TIndex _columnCount;

    std::span<const TIndex> _rowExtents;

    std::span<const TIndex> _columnIndices;

    std::span<const TMatrixValue> _entries;

    mutable OptionalRef<mp::ThreadPoolBase> _maybeThreads;
}; // class CSROperator


} // namespace cie::linalg
