#pragma once

// --- Linalg Includes ---
#include "packages/solvers/inc/LinearOperator.hpp"
#include "packages/solvers/inc/DefaultSpace.hpp"

// --- STL Includes ---
#include <memory>


namespace cie::linalg {


template <LinalgSpaceLike TSpace>
class JacobiOperator : public LinearOperator<TSpace> {
public:
    JacobiOperator()
    requires std::is_default_constructible_v<typename TSpace::Vector>;

    JacobiOperator(
        typename TSpace::Vector&& rInverseDiagonal,
        std::shared_ptr<const TSpace> pSpace) noexcept;

    void product(
        typename TSpace::ConstVectorView in,
        typename TSpace::Value scale,
        typename TSpace::VectorView out) const override;

private:
    std::shared_ptr<const TSpace> _pSpace;

    typename TSpace::Vector _inverseDiagonal;
}; // class JacobiOperator


template <
    class TValue,
    class TIndex,
    class TMatrixValue>
JacobiOperator<DefaultSpace<TValue,tags::SMP>> makeJacobiOperator(
    std::span<const TIndex> rowExtents,
    std::span<const TIndex> columnIndices,
    std::span<const TMatrixValue> entries,
    std::shared_ptr<const DefaultSpace<TValue,tags::SMP>> pSpace);


} // namespace cie::linalg
