// --- Linalg Includes ---
#include "packages/solvers/inc/JacobiOperator.hpp"
#include "packages/solvers/inc/DefaultSpace.hpp"
#include "packages/solvers/inc/SYCLSpace.hpp"

// --- Utility Includes ---
#include "packages/macros/inc/checks.hpp"

// --- STL Includes ---
#include <format>



namespace cie::linalg {


template <LinalgSpaceLike TS>
JacobiOperator<TS>::JacobiOperator(
    std::shared_ptr<TS> pSpace,
    std::size_t size,
    std::shared_ptr<LinearOperator<TS>> pLhs,
    std::shared_ptr<LinearOperator<TS>> pInverseDiagonal,
    std::size_t iterations,
    typename TS::Value relaxation)
        :   _pLhs(pLhs),
            _pInverseDiagonal(pInverseDiagonal),
            _iterations(iterations),
            _relaxation(relaxation),
            _pSpace(pSpace),
            _previous(pSpace->makeVector(size)),
            _residual(pSpace->makeVector(size)),
            _memory(pSpace->makeVector(size))
{}


template <LinalgSpaceLike TS>
void JacobiOperator<TS>::product(
    typename TS::Value inScale,
    typename TS::ConstVectorView in,
    typename TS::Value outScale,
    typename TS::VectorView out) {
        CIE_CHECK(
            _pSpace->size(in) == _pSpace->size(out) && _pSpace->size(out) == _pSpace->size(_previous),
            std::format(
                "inconsistent input vector sizes"))

        CIE_BEGIN_EXCEPTION_TRACING
            _pSpace->assign(_memory, out);
            _pSpace->fill(out, 0);
            //_pSpace->assign(_residual, in);

            for (std::size_t iIteration=0ul; iIteration<_iterations; ++iIteration) {
                _pSpace->assign(_previous, out);
                _pSpace->assign(_residual, in);
                _pLhs->product(1, out, -1, _residual);
                _pInverseDiagonal->product(1, _residual, _relaxation, out);
            } // for iIteration in range(_iterations)

            _pSpace->scale(out, outScale);
            _pSpace->add(out, _memory, inScale);
        CIE_END_EXCEPTION_TRACING
}


template class JacobiOperator<DefaultSpace<float>>;
template class JacobiOperator<DefaultSpace<double>>;


#ifdef CIE_ENABLE_SYCL

template class JacobiOperator<SYCLSpace<float>>;
template class JacobiOperator<SYCLSpace<double>>;

#endif


} // namespace cie::linalg
