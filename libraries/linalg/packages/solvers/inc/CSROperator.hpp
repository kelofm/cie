#pragma once

// --- Linalg Includes ---
#include "packages/solvers/inc/LinearOperator.hpp"
#include "packages/solvers/inc/DefaultSpace.hpp"
#include "packages/utilities/inc/CSRView.hpp"

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
        CSRView<const TMatrixValue,const TIndex> lhs,
        OptionalRef<mp::ThreadPoolBase> rMaybeThreads = {});

    /// @copydoc LinearOperator::product
    void product(
        typename Space::Value inScale,
        typename Space::ConstVectorView in,
        typename Space::Value outScale,
        typename Space::VectorView out) override;

protected:
    CSRView<const TMatrixValue,const TIndex> _lhs;

    OptionalRef<mp::ThreadPoolBase> _maybeThreads;
}; // class CSROperator


} // namespace cie::linalg
