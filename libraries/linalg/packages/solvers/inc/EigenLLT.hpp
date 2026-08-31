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
class EigenLLT : public LinearOperator<DefaultSpace<TValue>> {
private:
    using Space = DefaultSpace<TValue>;

public:
    EigenLLT();

    EigenLLT(EigenLLT&&) noexcept;

    EigenLLT(const EigenLLT&) = delete;

    EigenLLT& operator=(EigenLLT&&) noexcept;

    EigenLLT& operator=(const EigenLLT&) = delete;

    ~EigenLLT();

    EigenLLT(
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
}; // class EigenLLT


} // namespace cie::fem
