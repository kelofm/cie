#pragma once

// --- Linalg Includes ---
#include "packages/solvers/inc/LinearOperator.hpp"
#include "packages/solvers/inc/DefaultSpace.hpp"
#include "packages/utilities/inc/CSRView.hpp"

// --- STL Includes ---
#include <memory>


namespace cie::linalg {


template <LinalgSpaceLike TSpace>
class DiagonalOperator : public LinearOperator<TSpace> {
public:
    DiagonalOperator()
    requires std::is_default_constructible_v<typename TSpace::Vector>;

    DiagonalOperator(
        typename TSpace::Vector&& rInverseDiagonal,
        std::shared_ptr<const TSpace> pSpace) noexcept;

    void product(
        typename TSpace::Value inScale,
        typename TSpace::ConstVectorView in,
        typename TSpace::Value outScale,
        typename TSpace::VectorView out) override;

private:
    std::shared_ptr<const TSpace> _pSpace;

    typename TSpace::Vector _buffer, _inverseDiagonal;
}; // class DiagonalOperator


template <
    class TValue,
    class TIndex,
    class TMatrixValue>
DiagonalOperator<DefaultSpace<TValue,tags::SMP>> makeDiagonalOperator(
    CSRView<const TMatrixValue,const TIndex> matrix,
    std::shared_ptr<const DefaultSpace<TValue,tags::SMP>> pSpace);


} // namespace cie::linalg
