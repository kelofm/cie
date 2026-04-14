#pragma once

//--- Linalg Includes ---
#include "packages/solvers/inc/ConjugateGradients.hpp"

// --- Utility Includes ---
#include "packages/macros/inc/checks.hpp"
#include "packages/macros/inc/exceptions.hpp"
#include "packages/maths/inc/Comparison.hpp"
#include "packages/logging/inc/LoggerSingleton.hpp"
#include "packages/logging/inc/LogBlock.hpp"

// --- STL Includes ---
#include <format>
#include <cmath>
#include <optional>


namespace cie::linalg {


template <LinalgSpaceLike TS>
ConjugateGradients<TS>::ConjugateGradients()
requires std::is_default_constructible_v<TS>
    : ConjugateGradients(
        nullptr,
        std::make_shared<TS>(),
        nullptr,
        {},
        0)
{}


template <LinalgSpaceLike TS>
ConjugateGradients<TS>::ConjugateGradients(
    std::shared_ptr<LinearOperator<TS>> pLhs,
    std::shared_ptr<TS> pSpace,
    std::shared_ptr<LinearOperator<TS>> pPreconditioner,
    Statistics settings,
    int verbosity)
    :   IterativeSolver<TS>(settings),
        _pLhs(pLhs),
        _pPreconditioner(pPreconditioner),
        _pSpace(pSpace),
        _verbosity(verbosity)
{}


template <LinalgSpaceLike TS>
void ConjugateGradients<TS>::product(
    typename TS::Value inScale,
    typename TS::ConstVectorView in,
    typename TS::Value outScale,
    typename TS::VectorView out) {
        const std::size_t systemSize = _pSpace->size(out);

        CIE_CHECK(
            _pSpace->size(in) == systemSize,
            std::format(
                "incompatible vector sizes {} != {}",
                _pSpace->size(in), systemSize))

        const Statistics settings = this->getConfiguration();
        std::string report;
        Statistics stats {
            .iterationCount = 0,
            .absoluteResidual = std::numeric_limits<typename TS::Value>::max(),
            .relativeResidual = std::numeric_limits<typename TS::Value>::max()};

        std::optional<utils::LogBlock> maybeLogBlock;
        if (2 <= _verbosity)
            maybeLogBlock.emplace("ConjugateGradients", utils::LoggerSingleton::get());

        typename TS::Vector solution = _pSpace->makeVector(systemSize);
        _pSpace->assign(solution, out);

        if (settings.iterationCount) {
            CIE_BEGIN_EXCEPTION_TRACING
                // Define buffers.
                typename TS::Vector search = _pSpace->makeVector(systemSize);
                typename TS::Vector searchProduct = _pSpace->makeVector(systemSize);
                _pSpace->fill(searchProduct, 0);
                std::optional<typename TS::Vector> maybePreconditionedResidual;

                // Compute the initial residual.
                typename TS::Vector residual = _pSpace->makeVector(systemSize);
                _pSpace->assign(residual, in);
                _pLhs->product(1, solution, -1, residual);

                utils::Comparison<typename TS::Value> comparison;

                // Early exit if the initial residual satisfies the convergence criterion.
                typename TS::Value residualNorm = _pSpace->innerProduct(residual, residual);
                typename TS::Value preconditionedNorm = 0;
                const typename TS::Value initialResidualNorm = std::sqrt(residualNorm);

                stats.absoluteResidual = initialResidualNorm;
                stats.relativeResidual = 1;

                if (comparison.less(stats.absoluteResidual, settings.absoluteResidual) || comparison.less(stats.relativeResidual, settings.relativeResidual)) {
                    if (!this->makeIterationReport(
                        report,
                        _verbosity,
                        ConjugateGradients::ReportType::Termination,
                        stats,
                        settings).empty())
                            maybeLogBlock.value().log(report);

                    this->report(stats);
                    return;
                }

                // Compute the initial search direction.
                if (_pPreconditioner) {
                    maybePreconditionedResidual = _pSpace->makeVector(systemSize);
                    _pPreconditioner->product(0, residual, 1, *maybePreconditionedResidual);
                    preconditionedNorm = _pSpace->innerProduct(residual, *maybePreconditionedResidual);
                    _pSpace->assign(search, *maybePreconditionedResidual);
                    CIE_CHECK(
                        static_cast<typename TS::Value>(0) < preconditionedNorm,
                        std::format(
                            "preconditioned norm vanished in iteration {} (r^T @ z = 0)",
                            0))
                } else {
                    preconditionedNorm = residualNorm;
                    _pSpace->assign(search, residual);
                    CIE_CHECK(
                        static_cast<typename TS::Value>(0) < preconditionedNorm,
                        std::format(
                            "residual norm vanished in iteration {} (r^T @ r = 0)",
                            0))
                }

                for (; stats.iterationCount<settings.iterationCount; ++stats.iterationCount) {
                    // Compute part of the denominator of the search scale.
                    _pLhs->product(0, search, 1, searchProduct);

                    // Compute the search scale.
                    typename TS::Value searchScale = preconditionedNorm;
                    {
                        const typename TS::Value denominator = _pSpace->innerProduct(
                            search,
                            searchProduct);

                        CIE_CHECK(
                            static_cast<TS::Value>(0) < denominator,
                            std::format(
                                "p^T @ A @ p vanished in iteration {}",
                                stats.iterationCount))
                        searchScale /= denominator;
                    }

                    // Update the solution and residual.
                    _pSpace->add(solution, search, searchScale);
                    _pSpace->add(residual, searchProduct, -searchScale);

                    // Update the search direction.
                    residualNorm = _pSpace->innerProduct(residual, residual);
                    const typename TS::Value previousPreconditionedNorm = preconditionedNorm;

                    if (maybePreconditionedResidual) {
                        _pPreconditioner->product(0, residual, 1, *maybePreconditionedResidual);
                        preconditionedNorm = _pSpace->innerProduct(residual, *maybePreconditionedResidual);
                    } else {
                        preconditionedNorm = residualNorm;
                    }

                    searchScale = preconditionedNorm / previousPreconditionedNorm;
                    _pSpace->scale(search, searchScale);
                    _pSpace->add(search, maybePreconditionedResidual ? *maybePreconditionedResidual : residual, 1);

                    // Check whether the convergence criterion is satisfied.
                    stats.absoluteResidual = std::sqrt(residualNorm);
                    stats.relativeResidual = stats.absoluteResidual / initialResidualNorm;

                    if (!this->makeIterationReport(
                        report,
                        _verbosity,
                        ConjugateGradients::ReportType::Iteration,
                        stats,
                        settings).empty())
                            maybeLogBlock.value().log(report);

                    if (comparison.less(stats.absoluteResidual, settings.absoluteResidual) || comparison.less(stats.relativeResidual, settings.relativeResidual))
                        break;
                } // for stats.iterationCount in range(settings.iterationCount)
            CIE_END_EXCEPTION_TRACING
        }

        if (inScale == static_cast<typename TS::Value>(1)) {
            _pSpace->add(out, solution, outScale);
        } else if (inScale == static_cast<typename TS::Value>(0)) {
            _pSpace->assign(out, solution);
            _pSpace->scale(out, outScale);
        } else {
            _pSpace->scale(out, inScale);
            _pSpace->add(out, solution, outScale);
        }

        this->report(stats);
}


} // namespace cie::linalg
