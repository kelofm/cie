#pragma once

// --- Linalg Includes ---
#include "packages/constrained_solvers/inc/ConstrainedOperator.hpp"
#include "packages/solvers/inc/DefaultSpace.hpp"
#include "packages/solvers/inc/LinearSolverFactory.hpp"

// --- Utility Includes ---
#include "packages/io/inc/json.hpp"
#include "packages/stl_extension/inc/OptionalRef.hpp"
#include "packages/compile_time/packages/concepts/inc/basic_concepts.hpp"


namespace cie::linalg {


template <NumericLike TValue, IntegerLike TIndex>
class PenaltyConstrainedOperator
    : public ConstrainedOperator<DefaultSpace<TValue>,DefaultSpace<TIndex>> {
private:
    using Base = ConstrainedOperator<DefaultSpace<TValue>,DefaultSpace<TIndex>>;

public:
    using typename Base::ScalarSpace;

    using typename Base::IndexSpace;

    using typename Base::Scalar;

    using typename Base::Index;

    PenaltyConstrainedOperator() = default;

    PenaltyConstrainedOperator(PenaltyConstrainedOperator&&) noexcept;

    PenaltyConstrainedOperator& operator=(PenaltyConstrainedOperator&&) noexcept;

    PenaltyConstrainedOperator(
        Ref<const cie::io::JSONObject> rConfiguration,
        std::shared_ptr<ScalarSpace> pScalarSpace,
        std::shared_ptr<IndexSpace> pIndexSpace,
        CSRView<const Scalar,const Index> constraintGradients,
        typename ScalarSpace::ConstVectorView constraintGaps,
        CSRView<const Scalar,const Index> lhs,
        OptionalRef<LinearSolverFactory<ScalarSpace,IndexSpace>> rMaybeSolverFactory = {});

    ~PenaltyConstrainedOperator() override;

    /// @copydoc LinearOperator::product
    void product(
        typename ScalarSpace::Value inScale,
        typename ScalarSpace::ConstVectorView in,
        typename ScalarSpace::Value outScale,
        typename ScalarSpace::VectorView out) override;

private:
    PenaltyConstrainedOperator(const PenaltyConstrainedOperator&) = delete;

    PenaltyConstrainedOperator& operator=(const PenaltyConstrainedOperator&) = delete;

    struct Impl;
    std::unique_ptr<Impl> _pImpl;
}; // class PenaltyConstrainedOperator


} // namespace cie::linalg
