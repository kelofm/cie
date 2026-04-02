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

    ConjugateGradients()
    requires std::is_default_constructible_v<TSpace>;

    ConjugateGradients(
        std::shared_ptr<LinearOperator<TSpace>> pLhs,
        std::shared_ptr<TSpace> pSpace,
        std::shared_ptr<LinearOperator<TSpace>> pPreconditioner = nullptr,
        Statistics settings = {},
        int verbosity = 0);

    /// @copydoc LinearOperator::product
    void product(
        typename TSpace::Value inScale,
        typename TSpace::ConstVectorView in,
        typename TSpace::Value outScale,
        typename TSpace::VectorView out) override;

protected:
    std::shared_ptr<LinearOperator<TSpace>> _pLhs, _pPreconditioner;

    std::shared_ptr<TSpace> _pSpace;

    int _verbosity;
}; // class ConjugateGradients


} // namespace cie::fem
