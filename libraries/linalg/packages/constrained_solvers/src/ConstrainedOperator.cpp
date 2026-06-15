// --- Linalg Includes ---
#include "packages/constrained_solvers/inc/ConstrainedOperator.hpp"
#include "packages/solvers/inc/LinalgSpace.hpp"
#include "packages/solvers/inc/DefaultSpace.hpp"

// --- Utility Includes ---
#include "packages/macros/inc/exceptions.hpp"
#include "packages/types/inc/Color.hpp"

// --- STL Includes ---
#include <iostream>


namespace cie::linalg {


template <class T>
struct ConstrainedStatusStream<T>::Impl {
    typename ConstrainedStatusStream<T>::Status configuration;

    typename ConstrainedStatusStream<T>::Status status;
}; // struct ConstrainedStatusStream::Impl


template <class T>
ConstrainedStatusStream<T>::ConstrainedStatusStream()
    : ConstrainedStatusStream(
        /*constraintConfiguration=*/Status(),
        /*linearSystemConfiguration=*/Status(),
        /*rMabyeStream=*/{},
        /*scalarComparison=*/{})
{}


template <class T>
ConstrainedStatusStream<T>::ConstrainedStatusStream(
    Status constraintConfiguration,
    Status linearSystemConfiguration,
    OptionalRef<std::ostream> rMaybeStream,
    utils::Comparison<T> scalarComparison)
        :   StatusStream<T>(
                linearSystemConfiguration,
                rMaybeStream,
                scalarComparison),
            _pImpl(new Impl {
                .configuration = constraintConfiguration,
                .status = Status()})
{}


template <class T>
ConstrainedStatusStream<T>::ConstrainedStatusStream(ConstrainedStatusStream&&) noexcept = default;


template <class T>
ConstrainedStatusStream<T>::~ConstrainedStatusStream() = default;


template <class T>
ConstrainedStatusStream<T>&
ConstrainedStatusStream<T>::operator=(ConstrainedStatusStream&&) noexcept = default;


template <class T>
void ConstrainedStatusStream<T>::report(
    StatusReportType reportType,
    Verbosity verbosity) {
        CIE_BEGIN_EXCEPTION_TRACING
            Ref<std::ostream> rStream = this->stream();

            const bool constraintsAbsoluteConverged = this->lessEqual(
                    this->constraintStatus().absoluteResidual,
                    this->constraintConfiguration().absoluteResidual,
                    this->scalarComparison());
            const bool constraintsRelativeConverged = this->lessEqual(
                    this->constraintStatus().relativeResidual,
                    this->constraintConfiguration().relativeResidual,
                    this->scalarComparison());
            const bool constraintsConverged = constraintsAbsoluteConverged || constraintsRelativeConverged;

            const bool linearSystemAbsoluteConverged = this->lessEqual(
                    this->linearSystemStatus().absoluteResidual,
                    this->linearSystemConfiguration().absoluteResidual,
                    this->scalarComparison());
            const bool linearSystemRelativeConverged = this->lessEqual(
                    this->linearSystemStatus().relativeResidual,
                    this->linearSystemConfiguration().relativeResidual,
                    this->scalarComparison());
            const bool linearSystemConverged = linearSystemAbsoluteConverged || linearSystemRelativeConverged;

            if (Verbosity::Termination <= verbosity) {
                // Write header if necessary.
                if (!this->constraintStatus().iterationCount.has_value()) {
                    rStream
                        << "+ ----------------------------------------- + ----------------------------------------- +\n"
                        << "|               constraints                 |             linear system                 |\n"
                        << "+ --------- + ------------- + ------------- + --------- + ------------- + ------------- +\n"
                        << "| iteration | abs. residual | rel. residual | iteration | abs. residual | rel. residual |\n"
                        << "+ --------- + ------------- + ------------- + --------- + ------------- + ------------- +\n";
                } // if write header

                // Write status.
                if (reportType == StatusReportType::Termination
                    || constraintsConverged
                    || Verbosity::Iteration <= verbosity) {
                        rStream << std::format(
                            "| {}{:>9}{} | {}{:>13}{} | {}{:>13}{} | {}{:>9}{} | {}{:>13}{} | {}{:>13}{} |\n",
                            constraintsConverged ? RGBAColor::TUMGreen.ANSI() : RGBAColor::TUMOrange.ANSI(),
                            this->constraintStatus().iterationCount.has_value()
                                ? std::format("{}", this->constraintStatus().iterationCount.value())
                                : "",
                            "\033[0m",
                            constraintsAbsoluteConverged
                                ? RGBAColor::TUMGreen.ANSI()
                                : RGBAColor::TUMOrange.ANSI(),
                            this->constraintStatus().absoluteResidual.has_value()
                                ? std::format("{:.4E}", this->constraintStatus().absoluteResidual.value())
                                : "",
                            "\033[0m",
                            constraintsRelativeConverged
                                ? RGBAColor::TUMGreen.ANSI()
                                : RGBAColor::TUMOrange.ANSI(),
                            this->constraintStatus().relativeResidual.has_value()
                                ? std::format("{:.4E}", this->constraintStatus().relativeResidual.value())
                                : "",
                            "\033[0m",
                            linearSystemConverged ? RGBAColor::TUMGreen.ANSI() : RGBAColor::TUMOrange.ANSI(),
                            this->linearSystemStatus().iterationCount.has_value()
                                ? std::format("{}", this->linearSystemStatus().iterationCount.value())
                                : "",
                            "\033[0m",
                            linearSystemAbsoluteConverged
                                ? RGBAColor::TUMGreen.ANSI()
                                : RGBAColor::TUMOrange.ANSI(),
                            this->linearSystemStatus().absoluteResidual.has_value()
                                ? std::format("{:.4E}", this->linearSystemStatus().absoluteResidual.value())
                                : "",
                            "\033[0m",
                            linearSystemRelativeConverged
                                ? RGBAColor::TUMGreen.ANSI()
                                : RGBAColor::TUMOrange.ANSI(),
                            this->linearSystemStatus().relativeResidual.has_value()
                                ? std::format("{:.4E}", this->linearSystemStatus().relativeResidual.value())
                                : "",
                            "\033[0m");
                }

                // Write footer if necessary.
                if (constraintsConverged || reportType == StatusReportType::Termination)
                    rStream << "+ --------- + ------------- + ------------- + --------- + ------------- + ------------- +\n";
            } // Verbosity::Termination <= verbosity
        CIE_END_EXCEPTION_TRACING
}


template <class T>
Ref<const typename ConstrainedStatusStream<T>::Status>
ConstrainedStatusStream<T>::constraintConfiguration() const {
    return _pImpl->configuration;
}


template <class T>
Ref<const typename ConstrainedStatusStream<T>::Status>
ConstrainedStatusStream<T>::constraintStatus() const {
    return _pImpl->status;
}


template <class T>
void ConstrainedStatusStream<T>::configureConstraints(Status configuration) {
    _pImpl->configuration = configuration;
}


template <class T>
void ConstrainedStatusStream<T>::submitConstraintStatus(Status status) {
    _pImpl->status = status;
}


template <LinalgSpaceLike TSS, LinalgSpaceLike TIS>
ConstrainedOperator<TSS,TIS>::ConstrainedOperator(
    std::shared_ptr<ScalarSpace> pScalarSpace,
    std::shared_ptr<IndexSpace> pIndexSpace,
    CSRView<const Scalar,const Index> constraintGradients,
    typename ScalarSpace::ConstVectorView constraintGaps)
        :   _pScalarSpace(pScalarSpace),
            _pIndexSpace(pIndexSpace),
            _constraintGradients(constraintGradients),
            _constraintGaps(constraintGaps)
{}


#define CIE_INSTANTIATE_CONSTRAINED_OPERATOR(TS, T, I)  \
    template class ConstrainedOperator<TS<T>,TS<I>>;

#define CIE_INSTANTIATE_CONSTRAINED_BASES(TS, T)    \
    template class ConstrainedStatusStream<T>;      \
    CIE_INSTANTIATE_CONSTRAINED_OPERATOR(TS, T, int)

CIE_INSTANTIATE_CONSTRAINED_BASES(DefaultSpace, float)

CIE_INSTANTIATE_CONSTRAINED_BASES(DefaultSpace, double)

#undef CIE_INSTANTIATE_CONSTRAINED_BASES
#undef CIE_INSTANTIATE_CONSTRAINED_OPERATOR


} // namespace cie::linalg
