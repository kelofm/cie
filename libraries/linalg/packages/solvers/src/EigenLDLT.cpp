// --- External Includes ---
#include "Eigen/Sparse"
#include "Eigen/src/SparseCholesky/SimplicialCholesky.h"

// --- Linalg Includes ---
#include "packages/solvers/inc/EigenLDLT.hpp"

// --- Utility Includes ---
#include "packages/macros/inc/exceptions.hpp"

// --- STL Includes ---
#include <format>


namespace cie::linalg {


template <class T, class I>
struct EigenLDLT<T,I>::Impl {
    using EigenMatrix = Eigen::SparseMatrix<T,Eigen::RowMajor,I>;

    using EigenAdaptor = Eigen::Map<const EigenMatrix>;

    Eigen::SimplicialLDLT<EigenMatrix> solver;

    std::vector<T> buffer;

    OptionalRef<mp::ThreadPoolBase> maybeThreads;

    int verbosity;
}; // struct Impl


template <class T, class I>
EigenLDLT<T,I>::EigenLDLT(EigenLDLT&&) noexcept = default;


template <class T, class I>
EigenLDLT<T,I>& EigenLDLT<T,I>::operator=(EigenLDLT&&) noexcept = default;


template <class T, class I>
EigenLDLT<T,I>::~EigenLDLT() = default;


template <class T, class I>
EigenLDLT<T,I>::EigenLDLT(
    CSRView<const T,const I> lhs,
    OptionalRef<mp::ThreadPoolBase> rMaybeThreads)
        : _pImpl(new Impl) {
            typename Impl::EigenAdaptor lhsAdaptor(
                lhs.rowCount(),
                lhs.columnCount(),
                lhs.entries().size(),
                lhs.rowExtents().data(),
                lhs.columnIndices().data(),
                lhs.entries().data());

            CIE_BEGIN_EXCEPTION_TRACING
                _pImpl->solver.compute(lhsAdaptor);
            CIE_END_EXCEPTION_TRACING

            // Check whether the factorization was successful.
            switch (_pImpl->solver.info()) {
                case Eigen::ComputationInfo::NoConvergence:
                    CIE_THROW(Exception, "EigenLDLT failed to factorize the input matrix")
                    break;
                case Eigen::ComputationInfo::InvalidInput:
                    CIE_THROW(Exception, "Invalid input provided to EigenLDLT")
                    break;
                case Eigen::ComputationInfo::NumericalIssue:
                    CIE_THROW(Exception, "EigenLDLT ran into a numerical issue while factorizing the provided matrix")
                case Eigen::ComputationInfo::Success:
                    break;
            } // switch solver.info

            _pImpl->buffer.resize(lhs.rowCount(), T(0));
            _pImpl->maybeThreads = rMaybeThreads;
}


template <class T, class I>
void EigenLDLT<T,I>::product(
    typename Space::Value inScale,
    typename Space::ConstVectorView in,
    typename Space::Value outScale,
    typename Space::VectorView out) {
        CIE_BEGIN_EXCEPTION_TRACING
            using EigenSpan        = Eigen::Map<const Eigen::Matrix<T,Eigen::Dynamic,1>>;
            using EigenMutableSpan = Eigen::Map<Eigen::Matrix<T,Eigen::Dynamic,1>>;

            EigenSpan inAdaptor(in.data(), in.size(), 1);
            EigenMutableSpan bufferAdaptor(_pImpl->buffer.data(), _pImpl->buffer.size(), 1);
            EigenMutableSpan outAdaptor(out.data(), out.size(), 1);

            if (inScale == static_cast<T>(0)) {
                outAdaptor = _pImpl->solver.solve(inAdaptor);
                DefaultSpace<T>(_pImpl->maybeThreads).scale(out, outScale);
            } else {
                DefaultSpace<T>(_pImpl->maybeThreads).scale(out, inScale);
                bufferAdaptor = _pImpl->solver.solve(inAdaptor);
                DefaultSpace<T>(_pImpl->maybeThreads).add(
                    out,
                    _pImpl->buffer,
                    outScale);
            }
        CIE_END_EXCEPTION_TRACING
}


#define CIE_INSTANTIATE_EIGEN_LLT(T, I) template class EigenLDLT<T,I>;

CIE_INSTANTIATE_EIGEN_LLT(float, int)
CIE_INSTANTIATE_EIGEN_LLT(double, int)

#undef CIE_INSTANTIATE_EIGEN_LLT


} // namespace cie::linalg
