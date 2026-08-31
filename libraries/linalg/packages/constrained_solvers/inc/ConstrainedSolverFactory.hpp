#pragma once

// --- Linalg Includes ---
#include "packages/constrained_solvers/inc/ConstrainedOperator.hpp"
#include "packages/solvers/inc/LinearSolverFactory.hpp"
#include "packages/solvers/inc/LinalgSpace.hpp"
#include "packages/utilities/inc/CSRView.hpp"

// --- Utility Includes ---
#include "packages/registry/inc/FactoryRegistry.hpp"


namespace cie::linalg {


template <LinalgSpaceLike TScalarSpace, LinalgSpaceLike TIndexSpace>
class ConstrainedSolverFactory
    : public FactoryRegistry<
        ConstrainedOperator<TScalarSpace,TIndexSpace>,
        std::shared_ptr<TScalarSpace>,
        std::shared_ptr<TIndexSpace>,
        CSRView<
            const typename TScalarSpace::Value,
            const typename TIndexSpace::Value>,
        typename TScalarSpace::ConstVectorView,
        CSRView<
            const typename TScalarSpace::Value,
            const typename TIndexSpace::Value>,
        OptionalRef<LinearSolverFactory<TScalarSpace,TIndexSpace>>> {
private:
    using Base = FactoryRegistry<
        ConstrainedOperator<TScalarSpace,TIndexSpace>,
        std::shared_ptr<TScalarSpace>,
        std::shared_ptr<TIndexSpace>,
        CSRView<
            const typename TScalarSpace::Value,
            const typename TIndexSpace::Value>,
        typename TScalarSpace::ConstVectorView,
        CSRView<
            const typename TScalarSpace::Value,
            const typename TIndexSpace::Value>,
        OptionalRef<LinearSolverFactory<TScalarSpace,TIndexSpace>>>;

public:
    using ScalarSpace = TScalarSpace;

    using IndexSpace = TIndexSpace;

    using Scalar = typename ScalarSpace::Value;

    using Index = typename IndexSpace::Value;

    void load() override;
}; // class ConstrainedSolverFactory


template <class TScalarSpace, class TIndexSpace>
using ConstrainedSolverFactorySingleton = typename ConstrainedSolverFactory<
    TScalarSpace,
    TIndexSpace>::template Singleton<ConstrainedSolverFactory<
        TScalarSpace,
        TIndexSpace>>;


} // namespace cie::linalg
