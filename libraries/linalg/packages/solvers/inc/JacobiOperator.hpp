#pragma once

// --- Linalg Includes ---
#include "packages/solvers/inc/LinearOperator.hpp"
#include "packages/solvers/inc/DefaultSpace.hpp"

// --- Utility Includes ---
#include "packages/stl_extension/inc/OptionalRef.hpp"
#include "packages/concurrency/inc/ThreadPoolBase.hpp"


namespace cie::linalg {


/// @brief Linear operator representing an iteration of a scaled Jacobi relaxation.
/// @details Computes @f[
///             r += D^{-1} (c - (L + U) b)
///          @f]
///          where
///          - @f$A = L + D + U@f$ is a lower-diagonal-upper decomposition of the input matrix @f$A@f$,
///          - @f$D@f$ is a diagonal matrix,
///          - @f$L@f$ is a strictly lower triangular matrix,
///          - @f$U@f$ is a strictly upper triangular matrix,
///          - @f$c@f$ is the right-hand-side vector,
///          - @f$b@f$ is the input (solution) vector.
template <class TIndex, class TValue, class TMatrixValue = TValue>
class JacobiOperator
    : public LinearOperator<DefaultSpace<TValue,tags::SMP>> {
private:
    using Space = DefaultSpace<TValue,tags::SMP>;

public:
    constexpr JacobiOperator() noexcept = default;

    JacobiOperator(
        TIndex columnCount,
        std::span<const TIndex> rowExtents,
        std::span<const TIndex> columnIndices,
        std::span<const TMatrixValue> entries,
        std::size_t iterations,
        TValue relaxation,
        std::shared_ptr<Space> pSpace);

    /// @copydoc LinearOperator::product
    void product(
        typename Space::Value inScale,
        typename Space::ConstVectorView in,
        typename Space::Value outScale,
        typename Space::VectorView out) override;

protected:
    TIndex _columnCount;

    std::span<const TIndex> _rowExtents;

    std::span<const TIndex> _columnIndices;

    std::span<const TMatrixValue> _entries;

    typename Space::Vector _previous;

    std::size_t _iterations;

    TValue _relaxation;

    std::shared_ptr<Space> _pSpace;
}; // class JacobiOperator


} // namespace cie::linalg
