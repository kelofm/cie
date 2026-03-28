#pragma once

// --- Linalg Includes ---
#include "packages/solvers/inc/LinearOperator.hpp"
#include "packages/solvers/inc/LinalgSpace.hpp"

// --- Utility Includes ---
#include "packages/types/inc/types.hpp"


namespace cie::linalg {


template <LinalgSpaceLike TSpace>
class IterativeSolver {
public:
    using Operator = LinearOperator<TSpace>;

    struct Statistics {
        std::size_t iterationCount              = 0ul;
        typename TSpace::Value absoluteResidual = static_cast<typename TSpace::Value>(0);
        typename TSpace::Value relativeResidual = static_cast<typename TSpace::Value>(0);
    }; // struct Statistics

    virtual Statistics solve(
        Ref<const Operator> rLhs,
        typename TSpace::ConstVectorView rhs,
        typename TSpace::VectorView result,
        Statistics settings = {}) const = 0;
}; // class IterativeSolver


} // namespace cie::linalg
