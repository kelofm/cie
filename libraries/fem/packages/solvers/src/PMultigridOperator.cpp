#ifdef tmp

// --- FEM Includes ---
#include "packages/solvers/inc/PMultigridOperator.hpp"
#include "packages/macros/inc/exceptions.hpp"

// --- Linalg Includes ---
#include "packages/solvers/inc/CSROperator.hpp"

// --- STL Includes ---
#include <vector>

namespace cie::fem {


template <linalg::LinalgSpaceLike TSS, linalg::LinalgSpaceLike TIS>
struct PMultigridOperator<TSS,TIS>::Impl {
    std::vector<std::shared_ptr<linalg::LinearOperator<TSS>>> _gridOperators;

    std::shared_ptr<TSS> _pScalarSpace;

    std::shared_ptr<TIS> _pIndexSpace;

    typename TIS::Value _maxIterations;

    typename TSS::Value _absoluteResidual;

    typename TSS::Value _relativeResidual;
}; // struct Impl


template <linalg::LinalgSpaceLike TSS, linalg::LinalgSpaceLike TIS>
PMultigridOperator<TSS,TIS>::PMultigridOperator(
    std::shared_ptr<TSS> pScalarSpace,
    std::shared_ptr<TIS> pIndexSpace,
    linalg::CSRView<const typename TSS::Value,const typename TIS::Value> lhs,
    Ref<const Assembler> rAssembler,
    Ref<const cie::io::JSONObject> rConfiguration)
        : _pImpl(new Impl) {
            CIE_BEGIN_EXCEPTION_TRACING
                _pImpl->_pScalarSpace = std::move(pScalarSpace);
                _pImpl->_pIndexSpace  = std::move(pIndexSpace);

                //
            CIE_END_EXCEPTION_TRACING
}


template <linalg::LinalgSpaceLike TSS, linalg::LinalgSpaceLike TIS>
PMultigridOperator<TSS,TIS>::PMultigridOperator(PMultigridOperator&&) noexcept = default;


template <linalg::LinalgSpaceLike TSS, linalg::LinalgSpaceLike TIS>
PMultigridOperator<TSS,TIS>& PMultigridOperator<TSS,TIS>::operator=(PMultigridOperator<TSS,TIS>&&) noexcept = default;


template <linalg::LinalgSpaceLike TSS, linalg::LinalgSpaceLike TIS>
PMultigridOperator<TSS,TIS>::~PMultigridOperator() = default;


template <linalg::LinalgSpaceLike TSS, linalg::LinalgSpaceLike TIS>
void PMultigridOperator<TSS,TIS>::product(
    typename TSS::Value inScale,
    typename TSS::ConstVectorView in,
    typename TSS::Value outScale,
    typename TSS::VectorView out) {
        CIE_BEGIN_EXCEPTION_TRACING

        CIE_END_EXCEPTION_TRACING
}


} // namespace cie::fem

#endif