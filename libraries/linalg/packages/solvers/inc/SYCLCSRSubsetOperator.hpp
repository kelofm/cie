#ifdef CIE_ENABLE_SYCL

#pragma once

// --- Linalg Includes ---
#include "packages/solvers/inc/LinearOperator.hpp"
#include "packages/solvers/inc/SYCLSpace.hpp"
#include "packages/utilities/inc/CSRView.hpp"

// --- STL Includes ---
#include <memory>


namespace cie::linalg {


/// @brief Linear operator representing a scaled product of a vector with a subset of a matrix in CSR format.
template <class TIndex, class TValue>
class SYCLCSRSubsetOperator
    : public LinearOperator<SYCLSpace<TValue>> {
private:
    using Space = SYCLSpace<TValue>;

    using IndexSpace = SYCLSpace<TIndex>;

public:
    constexpr SYCLCSRSubsetOperator() noexcept = default;

    SYCLCSRSubsetOperator(
        CSRView<const TValue,const TIndex> lhs,
        typename IndexSpace::ConstVectorView subset,
        std::shared_ptr<Space> pSpace,
        std::shared_ptr<IndexSpace> pIndexSpace);

    /// @copydoc LinearOperator::product
    void product(
        typename Space::Value inScale,
        typename Space::ConstVectorView in,
        typename Space::Value outScale,
        typename Space::VectorView out) override;

protected:
    CSRView<const TValue,const TIndex> _lhs;

    typename IndexSpace::ConstVectorView _subset;

    std::shared_ptr<Space> _pSpace;

    std::shared_ptr<IndexSpace> _pIndexSpace;

    std::size_t _subGroupSize;

    std::size_t _groupSize;
}; // class SYCLCSRSubsetOperator


} // namespace cie::linalg

#endif
