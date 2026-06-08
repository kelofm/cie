// --- External Includes ---
#include "Eigen/Sparse"
#include "Eigen/IterativeLinearSolvers"

// --- Linalg Includes ---
#include "packages/solvers/inc/EigenCG.hpp"

// --- Utility Includes ---
#include "packages/macros/inc/exceptions.hpp"

// --- STL Includes ---
#include <variant>
#include <format>


namespace cie::linalg {


template <class T, class I>
struct EigenCG<T,I>::Impl {
    using EigenMatrix = Eigen::SparseMatrix<T,Eigen::RowMajor,I>;

    using EigenAdaptor = Eigen::Map<const EigenMatrix>;

    std::variant<
        Eigen::ConjugateGradient<
            EigenMatrix,
            Eigen::Lower | Eigen::Upper,
            Eigen::IdentityPreconditioner>,
        Eigen::ConjugateGradient<
            EigenMatrix,
            Eigen::Lower | Eigen::Upper,
            Eigen::DiagonalPreconditioner<T>>
    > solver;

    std::vector<T> buffer;

    OptionalRef<mp::ThreadPoolBase> maybeThreads;

    typename EigenCG<T,I>::Statistics settings;

    int verbosity;
}; // struct Impl


template <class T, class I>
EigenCG<T,I>::EigenCG(EigenCG&&) noexcept = default;


template <class T, class I>
EigenCG<T,I>& EigenCG<T,I>::operator=(EigenCG&&) noexcept = default;


template <class T, class I>
EigenCG<T,I>::~EigenCG() = default;


template <class T, class I>
EigenCG<T,I>::EigenCG(
    CSRView<const T,const I> lhs,
    I maxIterations,
    T relativeResidual,
    std::string_view preconditionerName,
    OptionalRef<mp::ThreadPoolBase> rMaybeThreads)
        : _pImpl(new Impl) {
            typename Impl::EigenAdaptor lhsAdaptor(
                lhs.rowCount(),
                lhs.columnCount(),
                lhs.entries().size(),
                lhs.rowExtents().data(),
                lhs.columnIndices().data(),
                lhs.entries().data());
            if (preconditionerName == "identity") {
                CIE_BEGIN_EXCEPTION_TRACING
                    _pImpl->solver.template emplace<Eigen::ConjugateGradient<
                        typename Impl::EigenMatrix,
                        Eigen::Lower | Eigen::Upper,
                        Eigen::IdentityPreconditioner>>();
                CIE_END_EXCEPTION_TRACING
            } else if (preconditionerName == "diagonal") {
                CIE_BEGIN_EXCEPTION_TRACING
                    _pImpl->solver.template emplace<Eigen::ConjugateGradient<
                        typename Impl::EigenMatrix,
                        Eigen::Lower | Eigen::Upper,
                        Eigen::DiagonalPreconditioner<T>>>();
                CIE_END_EXCEPTION_TRACING
            } else
                CIE_THROW(
                    Exception,
                    std::format(
                        "unsupported preconditioner \"{}\"",
                        preconditionerName))

            CIE_BEGIN_EXCEPTION_TRACING
                std::visit(
                    [lhsAdaptor, maxIterations, relativeResidual] (auto& rSolver) {
                        rSolver.setMaxIterations(maxIterations);
                        rSolver.setTolerance(relativeResidual);
                        rSolver.compute(lhsAdaptor);},
                    _pImpl->solver);
            CIE_END_EXCEPTION_TRACING

            _pImpl->buffer.resize(lhs.rowCount(), T(0));
            _pImpl->maybeThreads = rMaybeThreads;
            _pImpl->settings = typename EigenCG::Statistics {
                .iterationCount = static_cast<std::size_t>(maxIterations),
                .absoluteResidual = T(0),
                .relativeResidual = relativeResidual};
}


template <class T, class I>
void EigenCG<T,I>::product(
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
                std::visit(
                    [inAdaptor, &outAdaptor] (const auto& rSolver) {
                        outAdaptor = rSolver.solve(inAdaptor);},
                    _pImpl->solver);
                DefaultSpace<T>(_pImpl->maybeThreads).scale(out, outScale);
            } else {
                DefaultSpace<T>(_pImpl->maybeThreads).scale(out, inScale);
                std::visit(
                    [inAdaptor, &bufferAdaptor] (const auto& rSolver) {
                        bufferAdaptor = rSolver.solve(inAdaptor);},
                    _pImpl->solver);
                DefaultSpace<T>(_pImpl->maybeThreads).add(
                    out,
                    _pImpl->buffer,
                    outScale);
            }

            std::string report;
            typename EigenCG::Statistics stats;

            std::visit(
                [&] (const auto& rSolver) {
                    stats = typename EigenCG::Statistics {
                        .iterationCount   = static_cast<std::size_t>(rSolver.iterations()),
                        .absoluteResidual = std::numeric_limits<T>::quiet_NaN(),
                        .relativeResidual = rSolver.error()};},
                _pImpl->solver);
            this->makeIterationReport(
                report,
                _pImpl->verbosity,
                EigenCG::ReportType::Termination,
                stats,
                _pImpl->settings);
        CIE_END_EXCEPTION_TRACING
}


#define CIE_INSTANTIATE_EIGEN_CG(T, I) template class EigenCG<T,I>;

CIE_INSTANTIATE_EIGEN_CG(float, int)
CIE_INSTANTIATE_EIGEN_CG(double, int)

#undef CIE_INSTANTIATE_EIGEN_CG


} // namespace cie::linalg
