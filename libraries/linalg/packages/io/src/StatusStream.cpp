// --- Linalg Includes ---
#include "packages/io/inc/StatusStream.hpp"

// --- Utility Includes --
#include "packages/macros/inc/exceptions.hpp"
#include "packages/types/inc/Color.hpp"

// --- STL Includes ---
#include <iostream>


namespace cie::linalg {


template <class T>
struct StatusStream<T>::Impl {
    Ptr<std::ostream> pStream;

    typename StatusStream<T>::Status configuration;

    typename StatusStream<T>::Status status;

    utils::Comparison<T> scalarComparison;
}; // struct StatusStream::Impl


template <class T>
StatusStream<T>::StatusStream()
    : StatusStream(Status())
{}


template <class T>
StatusStream<T>::StatusStream(
    Status linearSystemConfiguration,
    OptionalRef<std::ostream> rMaybeStream,
    utils::Comparison<T> scalarComparison)
        : _pImpl(new Impl {
            .pStream = rMaybeStream.has_value() ? &rMaybeStream.value() : &std::cout,
            .configuration = linearSystemConfiguration,
            .status = Status(),
            .scalarComparison = scalarComparison})
{}


template <class T>
StatusStream<T>::StatusStream(StatusStream&&) noexcept = default;


template <class T>
StatusStream<T>::~StatusStream() = default;


template <class T>
StatusStream<T>& StatusStream<T>::operator=(StatusStream&&) noexcept = default;


template <class T>
void StatusStream<T>::report(
    StatusReportType reportType,
    Verbosity verbosity) {
        CIE_BEGIN_EXCEPTION_TRACING
            Ref<std::ostream> rStream = *_pImpl->pStream;
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
                if (!this->linearSystemStatus().iterationCount.has_value()) {
                    rStream
                        << "+ --------- + ------------- + ------------- +\n"
                        << "| iteration | abs. residual | rel. residual |\n"
                        << "+ --------- + ------------- + ------------- +\n";
                } // if write header

                // Write status.
                if (reportType == StatusReportType::Termination
                    || linearSystemConverged
                    || Verbosity::Iteration <= verbosity) {
                        rStream << std::format(
                            "| {}{:>9}{} | {}{:>13}{} | {}{:>13}{} |\n",
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
                if (linearSystemConverged || reportType == StatusReportType::Termination)
                    rStream << "+ --------- + ------------- + ------------- +\n";
            } // Verbosity::Termination <= verbosity
        CIE_END_EXCEPTION_TRACING
}


template <class T>
Ref<const typename StatusStream<T>::Status>
StatusStream<T>::linearSystemConfiguration() const {
    return _pImpl->configuration;
}


template <class T>
Ref<const typename StatusStream<T>::Status>
StatusStream<T>::linearSystemStatus() const {
    return _pImpl->status;
}


template <class T>
void StatusStream<T>::configureLinearSystem(Status configuration) {
    _pImpl->configuration = configuration;
}


template <class T>
void StatusStream<T>::submitLinearSystemStatus(Status status) {
    _pImpl->status = status;
}


template <class T>
Ref<const utils::Comparison<T>>
StatusStream<T>::scalarComparison() const {
    return _pImpl->scalarComparison;
}


template <class T>
bool StatusStream<T>::lessEqual(
    std::optional<T> maybeStatus,
    std::optional<T> maybeConfiguration,
    Ref<const utils::Comparison<T>> rScalarComparison) {
        if (!maybeConfiguration.has_value()) return true;
        if (maybeStatus.has_value()) {
            if (rScalarComparison.equal(maybeStatus.value(), maybeConfiguration.value()))
                return true;
            if (rScalarComparison.less(maybeStatus.value(), maybeConfiguration.value()))
                return true;
        }
        return false;
}


template <class T>
Ref<std::ostream> StatusStream<T>::stream() {
    return *_pImpl->pStream;
}


#define CIE_INSTANTIATE_STATUS_STREAM(T) template class StatusStream<T>;

CIE_INSTANTIATE_STATUS_STREAM(double)

CIE_INSTANTIATE_STATUS_STREAM(float)

#undef CIE_INSTANTIATE_STATUS_STREAM

} // namespace cie::linalg
