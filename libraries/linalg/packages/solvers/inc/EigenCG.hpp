#pragma once

// --- Linalg Includes ---
#include "packages/solvers/inc/IterativeSolver.hpp"
#include "packages/solvers/inc/DefaultSpace.hpp"
#include "packages/utilities/inc/CSRView.hpp"

// --- Utility Includes ---
#include "packages/stl_extension/inc/OptionalRef.hpp"

// --- STL Includes ---
#include <memory>
#include <string_view>


namespace cie::linalg {


template <class TValue, class TIndex>
class EigenCG : public IterativeSolver<DefaultSpace<TValue>> {
private:
    using Space = DefaultSpace<TValue>;

    using Base = IterativeSolver<Space>;

public:
    using typename Base::Status;

    EigenCG();

    EigenCG(EigenCG&&) noexcept;

    EigenCG(const EigenCG&) = delete;

    EigenCG& operator=(EigenCG&&) noexcept;

    EigenCG& operator=(const EigenCG&) = delete;

    ~EigenCG();

    EigenCG(
        CSRView<const TValue,const TIndex> lhs,
        TIndex maxIterations,
        TValue relativeResidual,
        std::string_view preconditionerName = "diagonal",
        Verbosity verbosity = Verbosity::Warnings,
        OptionalRef<mp::ThreadPoolBase> rMaybeThreads = {});

    /// @copydoc LinearOperator::product
    void product(
        typename Space::Value inScale,
        typename Space::ConstVectorView in,
        typename Space::Value outScale,
        typename Space::VectorView out) override;

protected:
    struct Impl;
    std::unique_ptr<Impl> _pImpl;
}; // class EigenCG


} // namespace cie::fem
