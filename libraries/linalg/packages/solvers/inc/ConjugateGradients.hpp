#pragma once

// --- Linalg Includes ---
#include "packages/solvers/inc/IterativeSolver.hpp"
#include "packages/solvers/inc/LinearOperator.hpp"

// --- STL Includes ---
#include <memory>
#include <type_traits>


namespace cie::linalg {


template <LinalgSpaceLike TSpace>
class ConjugateGradients : public IterativeSolver<TSpace> {
private:
    using Base = IterativeSolver<TSpace>;

public:
    using typename Base::Statistics;

    using Operator = LinearOperator<TSpace>;

    ConjugateGradients()
    requires std::is_default_constructible_v<TSpace>;

    ConjugateGradients(
        std::shared_ptr<TSpace> pSpace,
        Ptr<const Operator> pPreconditioner = nullptr,
        int verbosity = 0);

    Statistics solve(
        Ref<const Operator> rLhs,
        typename TSpace::ConstVectorView rhs,
        typename TSpace::VectorView result,
        Statistics settings = {}) const override;

protected:
    std::shared_ptr<TSpace> _pSpace;

    Ptr<const Operator> _pPreconditioner;

    int _verbosity;
}; // class ConjugateGradients


} // namespace cie::fem
