// --- Linalg Includes ---
#include "packages/solvers/inc/IterativeSolver.hpp"


namespace cie::linalg {


template <LinalgSpaceLike TS>
IterativeSolver<TS>::IterativeSolver(Statistics settings)
    :   _settings(settings),
        _results()
{}


template <LinalgSpaceLike TS>
constexpr typename IterativeSolver<TS>::Statistics
IterativeSolver<TS>::getConfiguration() const noexcept {
    return _settings;
}


template <LinalgSpaceLike TS>
constexpr void
IterativeSolver<TS>::configure(Statistics settings) noexcept {
    _settings = settings;
}


template <LinalgSpaceLike TS>
constexpr std::optional<typename IterativeSolver<TS>::Statistics>
IterativeSolver<TS>::getStats() const noexcept {
    return _results;
}


template <LinalgSpaceLike TS>
constexpr void
IterativeSolver<TS>::report(Statistics stats) noexcept {
    _results = stats;
}


} // namespace cie::linalg
