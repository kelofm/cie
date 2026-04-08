#pragma once

// --- Linalg Includes ---
#include "packages/solvers/inc/LinearOperator.hpp"
#include "packages/solvers/inc/SYCLSpace.hpp"
#include "packages/utilities/inc/CSRView.hpp"

// --- STL Includes ---
#include <memory>


namespace cie::linalg {


/// @brief Linear operator representing a scaled matrix-vector product in CSR format.
template <class TIndex, class TValue, class TMaskIndex>
class SYCLMaskedCSROperator
    : public LinearOperator<SYCLSpace<TValue>> {
private:
    using Space = SYCLSpace<TValue>;

    using MaskSpace = SYCLSpace<TMaskIndex>;

public:
    constexpr SYCLMaskedCSROperator() noexcept = default;

    SYCLMaskedCSROperator(
        CSRView<const TValue,const TIndex> lhs,
        typename MaskSpace::ConstVectorView mask,
        TMaskIndex threshold,
        std::shared_ptr<Space> pSpace,
        std::shared_ptr<MaskSpace> pMaskSpace);

    /// @copydoc LinearOperator::product
    void product(
        typename Space::Value inScale,
        typename Space::ConstVectorView in,
        typename Space::Value outScale,
        typename Space::VectorView out) override;

protected:
    CSRView<const TValue,const TIndex> _lhs;

    typename MaskSpace::ConstVectorView _mask;

    TMaskIndex _threshold;

    std::shared_ptr<Space> _pSpace;

    std::shared_ptr<MaskSpace> _pMaskSpace;

    std::size_t _subGroupSize;

    std::size_t _groupSize;
}; // class SYCLMaskedCSROperator


} // namespace cie::linalg
