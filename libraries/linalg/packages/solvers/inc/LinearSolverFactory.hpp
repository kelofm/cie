#pragma once

// --- Linalg Includes ---
#include "packages/solvers/inc/LinalgSpace.hpp"
#include "packages/solvers/inc/LoggedOperator.hpp"
#include "packages/utilities/inc/CSRView.hpp"

// --- Utility Includes ---
#include "packages/registry/inc/FactoryRegistry.hpp"


namespace cie::linalg {


template <LinalgSpaceLike TScalarSpace, LinalgSpaceLike TIndexSpace>
class LinearSolverFactory
    : public FactoryRegistry<
        LoggedOperator<TScalarSpace>,
        std::shared_ptr<TScalarSpace>,
        std::shared_ptr<TIndexSpace>,
        CSRView<
            const typename TScalarSpace::Value,
            const typename TIndexSpace::Value>> {
private:
    using Base = FactoryRegistry<
        LoggedOperator<TScalarSpace>,
        std::shared_ptr<TScalarSpace>,
        std::shared_ptr<TIndexSpace>,
        CSRView<
            const typename TScalarSpace::Value,
            const typename TIndexSpace::Value>>;

public:
    using ScalarSpace = TScalarSpace;

    using IndexSpace = TIndexSpace;

    using Scalar = typename ScalarSpace::Value;

    using Index = typename IndexSpace::Value;

    //using Singleton = typename Base::template Singleton<LinearSolverFactory>;

    void load() override;
}; // class LinearSolverFactory


template <class TScalarSpace, class TIndexSpace>
using LinearSolverFactorySingleton = typename LinearSolverFactory<
    TScalarSpace,
    TIndexSpace>::template Singleton<LinearSolverFactory<
        TScalarSpace,
        TIndexSpace>>;


} // namespace cie::linalg
