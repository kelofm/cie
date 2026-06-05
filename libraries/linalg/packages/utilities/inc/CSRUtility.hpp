#pragma once

// --- Utility Includes ---
#include "packages/concurrency/inc/ThreadPoolBase.hpp"
#include "packages/stl_extension/inc/OptionalRef.hpp"

// --- Linalg Includes ---
#include "packages/utilities/inc/CSRView.hpp"


namespace cie::linalg {


template <class TValue, class TIndex>
class CSRUtility {
public:
    CSRUtility() noexcept = default;

    CSRUtility(OptionalRef<mp::ThreadPoolBase> rMaybeThreads) noexcept;

    void sum(
        CSRView<TValue,const TIndex> left,
        CSRView<const TValue,const TIndex> right,
        TValue scale = static_cast<TValue>(1));

    void sum(
        Ref<std::vector<TIndex>> rTargetRowExtents,
        Ref<std::vector<TIndex>> rTargetColumnIndices,
        Ref<std::vector<TValue>> rTargetEntries,
        CSRView<const TValue,const TIndex> right,
        TValue targetScale = static_cast<TValue>(1),
        TValue rightScale = static_cast<TValue>(1));

    /// @brief Computes a symmetric matrix' product with itself.
    void symmetricProduct(
        Ref<std::vector<TIndex>> rTargetRowExtents,
        Ref<std::vector<TIndex>> rTargetColumnIndices,
        Ref<std::vector<TValue>> rTargetEntries,
        CSRView<const TValue,const TIndex> matrix,
        TValue scale = static_cast<TValue>(1));

    /// @brief Check whether two sorted arrays have at least one item in common.
    [[nodiscard]] bool hasSparseOverlap(
        std::span<const TIndex> left,
        std::span<const TIndex> right);

    /// @brief Compute the inner product of two sparse arrays.
    [[nodiscard]] TValue sparseInnerProduct(
        TIndex iRow,
        TIndex iColumn,
        Ref<const CSRView<const TValue, const TIndex>> rMatrix);

private:
    OptionalRef<mp::ThreadPoolBase> _rMaybeThreads;
}; // class CSRUtility


} // namespace cie::linalg
