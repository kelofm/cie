#pragma once

// --- Linalg Includes ---
#include "packages/solvers/inc/LinearOperator.hpp"
#include "packages/solvers/inc/SYCLSpace.hpp"
#include "packages/utilities/inc/CSRView.hpp"

// --- STL Includes ---
#include <memory>


namespace cie::linalg {


/// @brief Linear operator representing a scaled matrix-vector product in CSR format.
template <class TIndex, class TValue, class TMatrixValue = TValue>
class SYCLCSROperator
    : public LinearOperator<SYCLSpace<TValue>> {
private:
    using Space = SYCLSpace<TValue>;

public:
    constexpr SYCLCSROperator() noexcept = default;

    SYCLCSROperator(
        CSRView<const TMatrixValue,const TIndex> lhs,
        std::shared_ptr<SYCLSpace<TValue>> pSpace);

    /// @copydoc LinearOperator::product
    void product(
        typename Space::Value inScale,
        typename Space::ConstVectorView in,
        typename Space::Value outScale,
        typename Space::VectorView out) override;

protected:
    CSRView<const TMatrixValue,const TIndex> _lhs;

    std::shared_ptr<SYCLSpace<TValue>> _pSpace;

    std::size_t _subGroupSize;

    std::size_t _groupSize;
}; // class SYCLCSROperator


} // namespace cie::linalg
