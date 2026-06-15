#pragma once

// --- Linalg Includes ---
#include "packages/solvers/inc/LinearOperator.hpp"
#include "packages/solvers/inc/DefaultSpace.hpp"
#include "packages/utilities/inc/CSRView.hpp"

// --- Utility Includes ---
#include "packages/stl_extension/inc/OptionalRef.hpp"

// --- STL Includes ---
#include <memory>


namespace cie::linalg {


template <class TValue, class TIndex>
class EigenLDLT : public LinearOperator<DefaultSpace<TValue>> {
private:
    using Space = DefaultSpace<TValue>;

public:
    EigenLDLT();

    EigenLDLT(EigenLDLT&&) noexcept;

    EigenLDLT(const EigenLDLT&) = delete;

    EigenLDLT& operator=(EigenLDLT&&) noexcept;

    EigenLDLT& operator=(const EigenLDLT&) = delete;

    ~EigenLDLT();

    EigenLDLT(
        CSRView<const TValue,const TIndex> lhs,
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
}; // class EigenLDLT


} // namespace cie::fem
