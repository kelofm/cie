#pragma once

// --- Linalg Includes ---
#include "packages/solvers/inc/LoggedOperator.hpp"
#include "packages/solvers/inc/LinalgSpace.hpp"


namespace cie::linalg {


template <LinalgSpaceLike TSpace>
class IterativeSolver : public LoggedOperator<TSpace> {
public:
    using StreamLogger = StatusStream<typename TSpace::Value>;

    using Status = typename StreamLogger::Status;

    IterativeSolver() = default;

    explicit IterativeSolver(Ref<const Status> rConfiguration);

    [[nodiscard]] Status configuration() const;

    [[nodiscard]] Status status() const;

    void configure(Ref<const Status> rConfiguration);

protected:
    void updateStatus(Ref<const Status> rStatus);
}; // class IterativeSolver


} // namespace cie::linalg
