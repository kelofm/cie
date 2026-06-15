#pragma once

// --- Utility Includes ---
#include "packages/types/inc/types.hpp"
#include "packages/stl_extension/inc/OptionalRef.hpp"
#include "packages/maths/inc/Comparison.hpp"

// --- STL Includes ---
#include <memory>
#include <iosfwd>
#include <optional>


namespace cie::linalg {


enum class StatusReportType {
    Termination,
    Iteration
}; // enum class StatusReportType


enum class Verbosity {
    Silent      = 0,
    Warnings    = 1,
    Termination = 2,
    Iteration   = 3,
    Debug       = 4
}; // enum class Verbosity


template <class TScalar>
class StatusStream {
public:
    struct Status {
        std::optional<std::size_t>  iterationCount      = {};
        std::optional<TScalar>      absoluteResidual    = {};
        std::optional<TScalar>      relativeResidual    = {};
    }; // struct Status

    StatusStream();

    StatusStream(
        Status linearSystemConfiguration,
        OptionalRef<std::ostream> rMaybeStream = {},
        utils::Comparison<TScalar> scalarComparison = {});

    StatusStream(StatusStream&&) noexcept;

    virtual ~StatusStream();

    StatusStream& operator=(StatusStream&&) noexcept;

    virtual void report(
        StatusReportType reportType,
        Verbosity verbosity);

    [[nodiscard]] Ref<const Status> linearSystemConfiguration() const;

    [[nodiscard]] Ref<const Status> linearSystemStatus() const;

    void configureLinearSystem(Status configuration);

    void submitLinearSystemStatus(Status status);

    [[nodiscard]] Ref<const utils::Comparison<TScalar>> scalarComparison() const;

    [[nodiscard]] static bool lessEqual(
        std::optional<TScalar> maybeStatus,
        std::optional<TScalar> maybeConfiguration,
        Ref<const utils::Comparison<TScalar>> rScalarComparison);

    Ref<std::ostream> stream();

private:
    StatusStream(const StatusStream&) = delete;

    StatusStream& operator=(const StatusStream&) = delete;

    struct Impl;
    std::unique_ptr<Impl> _pImpl;
}; // class StatusStream


} // namespace cie::linalg
