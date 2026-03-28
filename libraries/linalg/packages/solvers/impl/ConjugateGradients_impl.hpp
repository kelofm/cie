#pragma once

//--- Linalg Includes ---
#include "packages/solvers/inc/ConjugateGradients.hpp"

// --- Utility Includes ---
#include "packages/macros/inc/checks.hpp"
#include "packages/macros/inc/exceptions.hpp"
#include "packages/maths/inc/Comparison.hpp"
#include "packages/logging/inc/LoggerSingleton.hpp"

// --- STL Includes ---
#include <format>
#include <cmath>
#include <optional>


namespace cie::linalg {


template <LinalgSpaceLike TS>
ConjugateGradients<TS>::ConjugateGradients()
requires std::is_default_constructible_v<TS>
    : ConjugateGradients(
        std::make_shared<TS>(),
        nullptr,
        0)
{}


template <LinalgSpaceLike TS>
ConjugateGradients<TS>::ConjugateGradients(
    std::shared_ptr<TS> pSpace,
        Ptr<const Operator> pPreconditioner,
    int verbosity)
    :
    _pSpace(pSpace),
    _pPreconditioner(pPreconditioner),
    _verbosity(verbosity)
{}


template <LinalgSpaceLike TS>
typename ConjugateGradients<TS>::Statistics
ConjugateGradients<TS>::solve(
    Ref<const Operator> rLhs,
    typename TS::ConstVectorView rhs,
    typename TS::VectorView result,
    Statistics settings) const {
        const std::size_t systemSize = _pSpace->size(result);

        CIE_CHECK(
            _pSpace->size(rhs) == systemSize,
            std::format(
                "incompatible vector sizes {} != {}",
                _pSpace->size(rhs), systemSize))

        Statistics output {
            .iterationCount = 0,
            .absoluteResidual = std::numeric_limits<typename TS::Value>::max(),
            .relativeResidual = std::numeric_limits<typename TS::Value>::max()};

        if (3 <= _verbosity) {
            utils::LoggerSingleton::get().log("ConjugateGradients");
            utils::LoggerSingleton::get().log("+ --------- + ----------------- + ----------------- +");
            utils::LoggerSingleton::get().log("| iteration | absolute residual | relative residual |");
            utils::LoggerSingleton::get().log("+ --------- + ----------------- + ----------------- +");
        }

        if (settings.iterationCount) {
            CIE_BEGIN_EXCEPTION_TRACING
                // Compute the initial residual.
                typename TS::Vector residual = _pSpace->makeVector(systemSize);
                _pSpace->assign(residual, rhs);
                rLhs.product(
                    result,
                    static_cast<typename TS::Value>(-1),
                    residual);

                utils::Comparison<typename TS::Value> comparison;

                // Early exit if the initial residual satisfies the convergence criterion.
                typename TS::Value residualNorm = _pSpace->innerProduct(residual, residual);
                typename TS::Value preconditionedNorm = 0;
                const typename TS::Value initialResidualNorm = std::sqrt(residualNorm);

                output.absoluteResidual = initialResidualNorm;
                output.relativeResidual = 1;
                if (comparison.less(output.absoluteResidual, settings.absoluteResidual) || comparison.less(output.relativeResidual, settings.relativeResidual))
                    return output;

                // Define buffers.
                typename TS::Vector search = _pSpace->makeVector(systemSize);
                typename TS::Vector searchProduct = _pSpace->makeVector(systemSize);
                std::optional<typename TS::Vector> maybePreconditionedResidual;

                // Compute the initial search direction.
                if (_pPreconditioner) {
                    maybePreconditionedResidual = _pSpace->makeVector(systemSize);
                    _pSpace->fill(*maybePreconditionedResidual, 0);
                    _pPreconditioner->product(residual, 1, *maybePreconditionedResidual);
                    preconditionedNorm = _pSpace->innerProduct(residual, *maybePreconditionedResidual);
                    _pSpace->assign(search, *maybePreconditionedResidual);
                } else {
                    preconditionedNorm = residualNorm;
                    _pSpace->assign(search, residual);
                }

                for (; output.iterationCount<settings.iterationCount; ++output.iterationCount) {
                    // Compute part of the denominator of the search scale.
                    _pSpace->fill(searchProduct, 0);
                    rLhs.product(search, 1, searchProduct);

                    // Compute the search scale.
                    typename TS::Value searchScale = preconditionedNorm / _pSpace->innerProduct(
                        search,
                        searchProduct);

                    // Update the solution and residual.
                    _pSpace->add(result, search, searchScale);
                    _pSpace->add(residual, searchProduct, -searchScale);

                    // Update the search direction.
                    residualNorm = _pSpace->innerProduct(residual, residual);
                    const typename TS::Value previousPreconditionedNorm = preconditionedNorm;

                    if (maybePreconditionedResidual) {
                        _pSpace->fill(*maybePreconditionedResidual, 0);
                        _pPreconditioner->product(residual, 1, *maybePreconditionedResidual);
                        preconditionedNorm = _pSpace->innerProduct(residual, *maybePreconditionedResidual);
                    } else {
                        preconditionedNorm = residualNorm;
                    }

                    searchScale = preconditionedNorm / previousPreconditionedNorm;
                    _pSpace->scale(search, searchScale);
                    _pSpace->add(search, maybePreconditionedResidual ? *maybePreconditionedResidual : residual, 1);

                    // Check whether the convergence criterion is satisfied.
                    output.absoluteResidual = std::sqrt(residualNorm);
                    output.relativeResidual = output.absoluteResidual / initialResidualNorm;

                    if (3 <= _verbosity)
                        utils::LoggerSingleton::get().log(std::format(
                            "| {:>9} | {:>17.5E} | {:>17.5E} |",
                            output.iterationCount, output.absoluteResidual, output.relativeResidual));

                    if (comparison.less(output.absoluteResidual, settings.absoluteResidual) || comparison.less(output.relativeResidual, settings.relativeResidual))
                        return output;
                } // for output.iterationCount in range(settings.iterationCount)
            CIE_END_EXCEPTION_TRACING
        }

        if (3 <= _verbosity)
            utils::LoggerSingleton::get().log("+ --------- + ----------------- + ----------------- +");

        return output;
}


} // namespace cie::linalg
