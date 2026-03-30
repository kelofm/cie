#pragma once

// --- Linalg Includes ---
#include "packages/solvers/inc/LinearOperator.hpp"
#include "packages/solvers/inc/DefaultSpace.hpp"

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
        typename TSpace::ConstVectorView in,
        typename TSpace::Value scale,
        typename TSpace::VectorView out) const override;

private:
    std::shared_ptr<const TSpace> _pSpace;

    typename TSpace::Vector _inverseDiagonal;
}; // class DiagonalOperator


template <
    class TValue,
    class TIndex,
    class TMatrixValue>
DiagonalOperator<DefaultSpace<TValue,tags::SMP>> makeDiagonalOperator(
    std::span<const TIndex> rowExtents,
    std::span<const TIndex> columnIndices,
    std::span<const TMatrixValue> entries,
    std::shared_ptr<const DefaultSpace<TValue,tags::SMP>> pSpace);


} // namespace cie::linalg
