// --- Linalg Includes ---
#include "packages/constrained_solvers/inc/ALMConstrainedOperator.hpp"
#include "packages/io/inc/StatusStream.hpp"
#include "packages/utilities/inc/CSRUtility.hpp"
#include "packages/solvers/inc/CSROperator.hpp"

// --- Utility Includes ---
#include "packages/macros/inc/exceptions.hpp"

// --- STL Includes ---
#include <cmath>


namespace cie::linalg {


template <NumericLike T, IntegerLike I>
struct ALMConstrainedOperator<T,I>::Impl {
    using Interface = ALMConstrainedOperator<T,I>;

    CSRView<const T,const I> unconstrainedLhs;

    typename Interface::IndexSpace::Vector constrainedRowExtents;

    typename Interface::IndexSpace::Vector constrainedColumnIndices;

    typename Interface::ScalarSpace::Vector constrainedEntries;

    std::shared_ptr<LoggedOperator<typename Interface::ScalarSpace>> pLinearOperator;

    std::shared_ptr<LinearOperator<typename Interface::ScalarSpace>> pConstraintGradientOperator;

    std::shared_ptr<LinearOperator<typename Interface::ScalarSpace>> pConstrainedLHSOperator;

    T penaltyFactor;

    std::shared_ptr<ConstrainedStatusStream<T>> pLogger;

    Verbosity verbosity;
}; // struct ALMConstrainedOperator::Impl


template <NumericLike T, IntegerLike I>
ALMConstrainedOperator<T,I>& ALMConstrainedOperator<T,I>::operator=(ALMConstrainedOperator&&) noexcept = default;


template <NumericLike T, IntegerLike I>
ALMConstrainedOperator<T,I>::~ALMConstrainedOperator() = default;


template <NumericLike T, IntegerLike I>
ALMConstrainedOperator<T,I>::ALMConstrainedOperator(
    Ref<const cie::io::JSONObject> rConfiguration,
    std::shared_ptr<ScalarSpace> pScalarSpace,
    std::shared_ptr<IndexSpace> pIndexSpace,
    CSRView<const Scalar,const Index> constraintGradients,
    typename ScalarSpace::ConstVectorView constraintGaps,
    CSRView<const Scalar,const Index> lhs,
    OptionalRef<LinearSolverFactory<ScalarSpace,IndexSpace>> rMaybeSolverFactory)
        :   Base(
                pScalarSpace,
                pIndexSpace,
                constraintGradients,
                constraintGaps),
            _pImpl(new Impl) {
                CIE_BEGIN_EXCEPTION_TRACING
                    // Parse configuration.
                    _pImpl->penaltyFactor = rConfiguration["penalty-factor"].as<double>();
                    _pImpl->verbosity = (Verbosity)rConfiguration["verbosity"].as<int>();
                CIE_END_EXCEPTION_TRACING
                CIE_BEGIN_EXCEPTION_TRACING
                    // Compute constrained LHS.
                    CSRUtility<T,I> utility(pScalarSpace->getThreads());
                    utility.symmetricProduct(
                        _pImpl->constrainedRowExtents,
                        _pImpl->constrainedColumnIndices,
                        _pImpl->constrainedEntries,
                        constraintGradients);
                    utility.sum(
                        _pImpl->constrainedRowExtents,
                        _pImpl->constrainedColumnIndices,
                        _pImpl->constrainedEntries,
                        lhs,
                        _pImpl->penaltyFactor);

                    // Construct an operator for the constrained LHS.
                    CSRView<const T,const I> constrainedLHS(
                            lhs.columnCount(),
                            _pImpl->constrainedRowExtents,
                            _pImpl->constrainedColumnIndices,
                            _pImpl->constrainedEntries);
                    Ref<LinearSolverFactory<ScalarSpace,IndexSpace>> rSolverFactory = rMaybeSolverFactory.has_value()
                        ? rMaybeSolverFactory.value()
                        : LinearSolverFactory<ScalarSpace,IndexSpace>::template Singleton<LinearSolverFactory<ScalarSpace,IndexSpace>>::get();
                    const auto maybeLinearOperator = rSolverFactory.make(
                        rConfiguration["linear-solver"],
                        pScalarSpace,
                        pIndexSpace,
                        constrainedLHS);
                    CIE_CHECK(
                        maybeLinearOperator.has_value(),
                        "linear solver could not be constructed")

                    _pImpl->pLinearOperator = maybeLinearOperator.value();

                    // Construct an operator for the constraint gradients.
                    _pImpl->pConstraintGradientOperator = std::make_shared<CSROperator<I,T,T>>(constraintGradients);
                    _pImpl->pConstrainedLHSOperator = std::make_shared<CSROperator<I,T,T>>(constrainedLHS);

                    // Set logger stream.
                    _pImpl->pLogger = std::make_shared<ConstrainedStatusStream<T>>(
                        typename StatusStream<T>::Status {
                            .iterationCount = rConfiguration["max-iterations"].as<std::size_t>(),
                            .absoluteResidual = rConfiguration["absolute-residual"].as<double>(),
                            .relativeResidual = rConfiguration["relative-residual"].as<double>()},
                        _pImpl->pLinearOperator->streamLogger().linearSystemConfiguration(),
                        _pImpl->pLinearOperator->streamLogger().stream(),
                        _pImpl->pLinearOperator->streamLogger().scalarComparison());
                    this->setLogger(_pImpl->pLogger);
                    _pImpl->pLinearOperator->setLogger(_pImpl->pLogger);
                CIE_END_EXCEPTION_TRACING
}


template <NumericLike T, IntegerLike I>
void ALMConstrainedOperator<T,I>::product(
    typename ScalarSpace::Value inScale,
    typename ScalarSpace::ConstVectorView in,
    typename ScalarSpace::Value outScale,
    typename ScalarSpace::VectorView out) {
        CIE_BEGIN_EXCEPTION_TRACING
            this->scalarSpace()->scale(out, inScale);

            // Allocate buffers.
            typename ScalarSpace::Vector constrainedRHS = this->scalarSpace()->makeVector(in.size());
            typename ScalarSpace::Vector solutionUpdate = this->scalarSpace()->makeVector(out.size());
            typename ScalarSpace::Vector solution       = this->scalarSpace()->makeVector(out.size());
            this->scalarSpace()->fill(solutionUpdate, 0);
            this->scalarSpace()->assign(solution, out);

            // Compute the initial constraint violation.
            typename ScalarSpace::Vector constraintResidual = this->scalarSpace()->makeVector(this->constraintGaps().size());
            this->scalarSpace()->assign(
                constraintResidual,
                this->constraintGaps());
            _pImpl->pConstraintGradientOperator->product(1, solution, 1, constraintResidual);
            const T initialConstraintResidualNorm = std::sqrt(this->scalarSpace()->innerProduct(
                constraintResidual,
                constraintResidual));

            // Compute the initial RHS.
            this->scalarSpace()->assign(
                constrainedRHS,
                this->constraintGaps());

            // Check for early exit possibilities.
            const typename StatusStream<T>::Status configuration = _pImpl->pLogger->constraintConfiguration();
            typename StatusStream<T>::Status status {
                .iterationCount = {},
                .absoluteResidual = initialConstraintResidualNorm,
                .relativeResidual = T(1)};

            _pImpl->pLogger->submitConstraintStatus(status);
            if (_pImpl->pLogger->lessEqual(
                    status.absoluteResidual,
                    configuration.absoluteResidual,
                    _pImpl->pLogger->scalarComparison())) {
                        _pImpl->pLogger->report(
                            StatusReportType::Termination,
                            _pImpl->verbosity);
                        return;
            } else {
                _pImpl->pLogger->report(
                    StatusReportType::Iteration,
                    _pImpl->verbosity);
            }
            CIE_CHECK(initialConstraintResidualNorm, "")

            // Run the constraint loop.
            for (std::size_t iConstraintIteration=0ul; iConstraintIteration<configuration.iterationCount.value(); ++iConstraintIteration) {
                // Update the RHS.
                _pImpl->pConstraintGradientOperator->product(
                    1, constraintResidual,
                    -_pImpl->penaltyFactor, constrainedRHS);

                // Solve the current linear system and apply the update.
                CIE_BEGIN_EXCEPTION_TRACING
                    _pImpl->pLinearOperator->product(0, constrainedRHS, 1, solutionUpdate);
                    this->scalarSpace()->add(solution, solutionUpdate, outScale);
                    _pImpl->pConstrainedLHSOperator->product(1, solutionUpdate, -1, constrainedRHS);
                CIE_END_EXCEPTION_TRACING

                // Compute the constraint violation.
                this->scalarSpace()->assign(
                    constraintResidual,
                    this->constraintGaps());
                _pImpl->pConstraintGradientOperator->product(1, solution, 1, constraintResidual);
                const T constraintResidualNorm = std::sqrt(this->scalarSpace()->innerProduct(
                    constraintResidual,
                    constraintResidual));

                // Check for convergence.
                status.iterationCount = iConstraintIteration;
                status.absoluteResidual = constraintResidualNorm;
                status.relativeResidual = constraintResidualNorm / initialConstraintResidualNorm;
                _pImpl->pLogger->submitConstraintStatus(status);
                if (iConstraintIteration == configuration.iterationCount.value() - 1 || _pImpl->pLogger->lessEqual(
                        status.absoluteResidual,
                        configuration.absoluteResidual,
                        _pImpl->pLogger->scalarComparison())) {
                            _pImpl->pLogger->report(
                                StatusReportType::Termination,
                                _pImpl->verbosity);
                            break;
                } else
                    _pImpl->pLogger->report(
                        StatusReportType::Iteration,
                        _pImpl->verbosity);
            } // for iConstraintIteration

            this->scalarSpace()->add(
                out,
                solution, outScale);
        CIE_END_EXCEPTION_TRACING
}


#define CIE_INSTANTIATE_PENALTY_OPERATOR(T, I)  \
    template class ALMConstrainedOperator<T,I>;

CIE_INSTANTIATE_PENALTY_OPERATOR(float, int)

CIE_INSTANTIATE_PENALTY_OPERATOR(double, int)

#undef CIE_INSTANTIATE_PENALTY_OPERATOR


} // namespace cie::linalg
