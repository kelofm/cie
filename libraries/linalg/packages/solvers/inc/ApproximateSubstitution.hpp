#pragma once

// --- Linalg Includes ---
#include "packages/solvers/inc/DefaultSpace.hpp"
#include "packages/solvers/inc/LinearOperator.hpp"
#include "packages/utilities/inc/CSRView.hpp"

// --- STL Includes ---
#include <memory>


namespace cie::linalg {


template <class TValue, class TIndex, class TMatrixValue = TValue>
class ApproximateSubstitution : public LinearOperator<DefaultSpace<TValue>> {
private:
    using Space = DefaultSpace<TValue>;

public:
    constexpr ApproximateSubstitution() noexcept = default;

    ApproximateSubstitution(
        CSRView<const TMatrixValue,const TIndex> lhs,
        std::size_t iterations,
        std::shared_ptr<Space> pSpace);

    void product(
        typename Space::Value inScale,
        typename Space::ConstVectorView in,
        typename Space::Value outScale,
        typename Space::VectorView out) override;

private:
    CSRView<const TMatrixValue,const TIndex> _lhs;

    std::shared_ptr<Space> _pSpace;

    std::size_t _iterations;

    typename Space::Vector _previous;
}; // class ApproximateSubstitution


} // namespace cie::linalg
