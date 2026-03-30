#pragma once

// --- Linalg Includes ---
#include "packages/solvers/inc/LinearOperator.hpp"
#include "packages/solvers/inc/LinalgSpace.hpp"

// --- STL Includes ---
#include <optional>


namespace cie::linalg {


template <LinalgSpaceLike TSpace>
class IterativeSolver : public LinearOperator<TSpace> {
public:
    struct Statistics {
        std::size_t iterationCount              = 0ul;
        typename TSpace::Value absoluteResidual = static_cast<typename TSpace::Value>(0);
        typename TSpace::Value relativeResidual = static_cast<typename TSpace::Value>(0);
    }; // struct Statistics

    constexpr IterativeSolver() noexcept = default;

    IterativeSolver(Statistics settings);

    [[nodiscard]] constexpr Statistics getConfiguration() const noexcept;

    constexpr void configure(Statistics settings) noexcept;

    [[nodiscard]] constexpr std::optional<Statistics> getStats() const noexcept;

protected:
    constexpr void report(Statistics stats) noexcept;

private:
    Statistics _settings;

    std::optional<Statistics> _results;
}; // class IterativeSolver


} // namespace cie::linalg

#include "packages/solvers/impl/IterativeSolver_impl.hpp"
