#pragma once

// --- FEM Includes ---
#include "packages/graph/inc/Assembler.hpp"

// --- LinAlg Includes ---
#include "packages/solvers/inc/LinearOperator.hpp"
#include "packages/solvers/inc/DefaultSpace.hpp"
#include "packages/solvers/inc/SYCLSpace.hpp"
#include "packages/utilities/inc/CSRView.hpp"

// --- Utility Includes ---
#include "packages/io/inc/json.hpp"

// --- STL Includes ---
#include <memory>


namespace cie::fem {


template <linalg::LinalgSpaceLike TScalarSpace, linalg::LinalgSpaceLike TIndexSpace>
class PMultigridOperator : public linalg::LinearOperator<TScalarSpace> {
public:
    PMultigridOperator() noexcept = default;

    PMultigridOperator(
        std::shared_ptr<TScalarSpace> pScalarSpace,
        std::shared_ptr<TIndexSpace> pIndexSpace,
        linalg::CSRView<const typename TScalarSpace::Value,const typename TIndexSpace::Value> lhs,
        Ref<const Assembler> rAssembler,
        Ref<const cie::io::JSONObject> rConfiguration);

    PMultigridOperator(PMultigridOperator&&) noexcept;

    PMultigridOperator(const PMultigridOperator&) = default;

    PMultigridOperator& operator=(PMultigridOperator&&) noexcept;

    PMultigridOperator& operator=(const PMultigridOperator&) = default;

    ~PMultigridOperator();

    /// @copydoc linalg::LinearOperator::product
    void product(
        typename TScalarSpace::Value inScale,
        typename TScalarSpace::ConstVectorView in,
        typename TScalarSpace::Value outScale,
        typename TScalarSpace::VectorView out) override;

private:
    struct Impl;
    std::unique_ptr<Impl> _pImpl;
}; // class PMultigridOperator


} // namespace cie::fem
