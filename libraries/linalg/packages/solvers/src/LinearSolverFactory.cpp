// --- Linalg Includes ---
#include "packages/solvers/inc/LinearSolverFactory.hpp"
#include "packages/solvers/inc/DefaultSpace.hpp"
#include "packages/solvers/inc/SYCLSpace.hpp"
#include "packages/solvers/inc/JacobiOperator.hpp"
#include "packages/solvers/inc/CSROperator.hpp"
#include "packages/solvers/inc/DiagonalOperator.hpp"
#include "packages/solvers/inc/ConjugateGradients.hpp"
#include "packages/solvers/inc/EigenCG.hpp"
#include "packages/solvers/inc/EigenLLT.hpp"
#include "packages/solvers/inc/EigenLDLT.hpp"
#include "packages/solvers/inc/SYCLCSROperator.hpp"

// --- Utility Includes ---
#include "packages/macros/inc/exceptions.hpp"
#include "packages/io/inc/json.hpp"
#include "packages/registry/impl/FactoryRegistry_impl.hpp"


namespace cie::linalg {


namespace detail {
template <class TSourceSpace, class TTargetSpace>
concept CrossAssignableSpaces
=  LinalgSpaceLike<TSourceSpace>
&& LinalgSpaceLike<TTargetSpace>
&& requires (
    typename TSourceSpace::ConstVectorView source,
    typename TTargetSpace::VectorView target,
    TTargetSpace& targetSpace) {
        {targetSpace.assign(target, source)} -> std::same_as<void>;
}; // concept CrossAssignableSpaces
} // namespace detail


template <
    LinalgSpaceLike TSourceScalarSpace,
    LinalgSpaceLike TSourceIndexSpace,
    LinalgSpaceLike TTargetScalarSpace,
    LinalgSpaceLike TTargetIndexSpace>
class CastSolver final : public LoggedOperator<TSourceScalarSpace> {
public:
    CastSolver() noexcept = default;

    CastSolver(
        CSRView<const typename TSourceScalarSpace::Value,const typename TSourceIndexSpace::Value> lhs,
        std::shared_ptr<TSourceScalarSpace> pSourceScalarSpace,
        std::shared_ptr<TSourceIndexSpace> pSourceIndexSpace,
        std::shared_ptr<TTargetScalarSpace> pTargetScalarSpace,
        std::shared_ptr<TTargetIndexSpace> pTargetIndexSpace,
        Ref<const cie::io::JSONObject> rConfiguration,
        Ref<const LinearSolverFactory<TTargetScalarSpace,TTargetIndexSpace>> rFactory)
            :   _pOperator(),
                _lhsRowExtents(pTargetIndexSpace->makeVector(lhs.rowExtents().size())),
                _lhsColumnIndices(pTargetIndexSpace->makeVector(lhs.columnIndices().size())),
                _lhsEntries(pTargetScalarSpace->makeVector(lhs.entries().size())),
                _pSourceScalarSpace(pSourceScalarSpace),
                _pTargetScalarSpace(pTargetScalarSpace) {
                    CIE_BEGIN_EXCEPTION_TRACING
                        this->copy<TSourceIndexSpace,TTargetIndexSpace>(
                            *pSourceIndexSpace,
                            *pTargetIndexSpace,
                            pSourceIndexSpace->view(lhs.rowExtents()),
                            pTargetIndexSpace->view(_lhsRowExtents));
                    CIE_END_EXCEPTION_TRACING
                    CIE_BEGIN_EXCEPTION_TRACING
                        this->copy<TSourceIndexSpace,TTargetIndexSpace>(
                            *pSourceIndexSpace,
                            *pTargetIndexSpace,
                            pSourceIndexSpace->view(lhs.columnIndices()),
                            pTargetIndexSpace->view(_lhsColumnIndices));
                    CIE_END_EXCEPTION_TRACING
                    CIE_BEGIN_EXCEPTION_TRACING
                        this->copy<TSourceScalarSpace,TTargetScalarSpace>(
                            *pSourceScalarSpace,
                            *pTargetScalarSpace,
                            pSourceScalarSpace->view(lhs.entries()),
                            pTargetScalarSpace->view(_lhsEntries));
                    CIE_END_EXCEPTION_TRACING
                    CIE_BEGIN_EXCEPTION_TRACING
                        CSRView<const typename TTargetScalarSpace::Value,const typename TTargetIndexSpace::Value> targetLhs(
                            lhs.columnCount(),
                            std::span<const typename TTargetIndexSpace::Value>(
                                pTargetIndexSpace->data(_lhsRowExtents),
                                pTargetIndexSpace->size(_lhsRowExtents)),
                            std::span<const typename TTargetIndexSpace::Value>(
                                pTargetIndexSpace->data(_lhsColumnIndices),
                                pTargetIndexSpace->size(_lhsColumnIndices)),
                            std::span<const typename TTargetScalarSpace::Value>(
                                pTargetScalarSpace->data(_lhsEntries),
                                pTargetScalarSpace->size(_lhsEntries)));
                        auto pMaybeOperator = rFactory.make(
                            rConfiguration,
                            pTargetScalarSpace,
                            pTargetIndexSpace,
                            targetLhs);
                        _pOperator = pMaybeOperator.value();
                        LoggedOperator<TSourceScalarSpace>::setLogger(_pOperator->streamLogger());
                    CIE_END_EXCEPTION_TRACING
    }


    void product(
        typename TSourceScalarSpace::Value inScale,
        typename TSourceScalarSpace::ConstVectorView in,
        typename TSourceScalarSpace::Value outScale,
        typename TSourceScalarSpace::VectorView out) override {
            CIE_BEGIN_EXCEPTION_TRACING
                typename TTargetScalarSpace::Vector targetIn = _pTargetScalarSpace->makeVector(_pSourceScalarSpace->size(in));
                this->copy<TSourceScalarSpace,TTargetScalarSpace>(
                    *_pSourceScalarSpace,
                    *_pTargetScalarSpace,
                    in,
                    targetIn);
                typename TTargetScalarSpace::Vector targetOut = _pTargetScalarSpace->makeVector(_pSourceScalarSpace->size(out));
                this->copy<TSourceScalarSpace,TTargetScalarSpace>(
                    *_pSourceScalarSpace,
                    *_pTargetScalarSpace,
                    out,
                    targetOut);
                _pOperator->product(
                    inScale,
                    targetIn,
                    outScale,
                    targetOut);
                this->copy<TTargetScalarSpace,TSourceScalarSpace>(
                    *_pTargetScalarSpace,
                    *_pSourceScalarSpace,
                    targetOut,
                    out);
            CIE_END_EXCEPTION_TRACING
    }

    void setLogger(Ref<const std::shared_ptr<StatusStream<typename TSourceScalarSpace::Value>>> rpLogger) override {
        LoggedOperator<TSourceScalarSpace>::setLogger(_pOperator->streamLogger());
        _pOperator->setLogger(rpLogger);
    }

private:
    /// @brief Copies an array from one space to another.
    /// @details Checks whether @p TTargetSpace can directly assign an array
    ///          from @p TSourceSpace. If it can, it calls the assignment.
    ///          If not, detours are taken through @ref DefaultSpace.
    template <class TSourceSpace, class TTargetSpace>
    void copy(
        TSourceSpace& rSourceSpace,
        TTargetSpace& rTargetSpace,
        typename TSourceSpace::ConstVectorView source,
        typename TTargetSpace::VectorView target) {
            CIE_CHECK(
                source.size() == target.size(),
                "incomatible array sizes " << source.size() << " != " << target.size())
            CIE_BEGIN_EXCEPTION_TRACING
                if constexpr (detail::CrossAssignableSpaces<TSourceSpace,TTargetSpace>) {
                    rTargetSpace.assign(target, source);
                } else {
                    DefaultSpace<typename TSourceSpace::Value> defaultSourceSpace;
                    auto buffer = defaultSourceSpace.makeVector(source.size());
                    rSourceSpace.assign(buffer, std::span<const typename TSourceSpace::Value>(source));
                    if constexpr (std::is_same_v<typename TSourceSpace::Value,typename TTargetSpace::Value>) {
                        rTargetSpace.assign(target, std::span<const typename TTargetSpace::Value>(buffer));
                    } else {
                        DefaultSpace<typename TTargetSpace::Value> defaultTargetSpace;
                        auto targetBuffer = defaultTargetSpace.makeVector(buffer.size());
                        std::copy(
                            buffer.begin(),
                            buffer.begin() + buffer.size(),
                            targetBuffer.begin());
                        rTargetSpace.assign(target, std::span<const typename TTargetSpace::Value>(targetBuffer));
                    }
                }
            CIE_END_EXCEPTION_TRACING
    }

    std::shared_ptr<LoggedOperator<TTargetScalarSpace>> _pOperator;

    typename TTargetIndexSpace::Vector _lhsRowExtents;

    typename TTargetIndexSpace::Vector _lhsColumnIndices;

    typename TTargetScalarSpace::Vector _lhsEntries;

    std::shared_ptr<TSourceScalarSpace> _pSourceScalarSpace;

    std::shared_ptr<TTargetScalarSpace> _pTargetScalarSpace;
}; // class CastSolver


template <class T, class I>
void registerInverseDiagonal(Ref<LinearSolverFactory<
        DefaultSpace<T>,
        DefaultSpace<I>>
    > rFactory) {
        CIE_BEGIN_EXCEPTION_TRACING
            cie::io::JSONObject schema(std::string {R"({"$ref" : "/cie/linalg/inverse-diagonal-operator"})"});
            const auto factory = [] (
                Ref<const cie::io::JSONObject>,
                std::shared_ptr<DefaultSpace<T>> pScalarSpace,
                std::shared_ptr<DefaultSpace<I>>,
                CSRView<const T,const I> lhs) -> std::shared_ptr<LoggedOperator<DefaultSpace<T>>> {
                    CIE_BEGIN_EXCEPTION_TRACING
                        return makeLoggedOperator<DiagonalOperator<DefaultSpace<T>>>(
                            makeDiagonalOperator<T,I,T>(
                                lhs,
                                pScalarSpace));
                    CIE_END_EXCEPTION_TRACING
                };
            rFactory.insert("inverse-diagonal-operator", schema, factory);
        CIE_END_EXCEPTION_TRACING
}


template <class T, class I>
void registerJacobi(Ref<LinearSolverFactory<
        DefaultSpace<T>,
        DefaultSpace<I>>
    > rFactory) {
        CIE_BEGIN_EXCEPTION_TRACING
            cie::io::JSONObject schema(std::string {R"({"$ref" : "/cie/linalg/jacobi"})"});
            const auto factory = [] (
                Ref<const cie::io::JSONObject> rConfiguration,
                std::shared_ptr<DefaultSpace<T>> pScalarSpace,
                std::shared_ptr<DefaultSpace<I>>,
                CSRView<const T,const I> lhs) -> std::shared_ptr<LoggedOperator<DefaultSpace<T>>> {
                    CIE_BEGIN_EXCEPTION_TRACING
                        auto pLHSOperator = std::make_shared<CSROperator<I,T,T>>(
                            lhs,
                            pScalarSpace->getThreads());
                        auto pDiagonalOperator = std::make_shared<DiagonalOperator<DefaultSpace<T>>>(
                            makeDiagonalOperator<T,I,T>(
                                lhs,
                                pScalarSpace));
                        return makeLoggedOperator<JacobiOperator<DefaultSpace<T>>>(
                            pScalarSpace,
                            lhs.rowCount(),
                            pLHSOperator,
                            pDiagonalOperator,
                            rConfiguration["iterations"].as<double>(),
                            rConfiguration["relaxation"].as<double>());
                    CIE_END_EXCEPTION_TRACING
                };
            rFactory.insert("jacobi", schema, factory);
        CIE_END_EXCEPTION_TRACING
}


template <class T, class I>
void registerCG(Ref<LinearSolverFactory<
        DefaultSpace<T>,
        DefaultSpace<I>>
    > rFactory) {
        CIE_BEGIN_EXCEPTION_TRACING
            cie::io::JSONObject schema(std::string {R"({"$ref" : "/cie/linalg/cg"})"});
            const auto factory = [&rFactory] (
                Ref<const cie::io::JSONObject> rConfiguration,
                std::shared_ptr<DefaultSpace<T>> pScalarSpace,
                std::shared_ptr<DefaultSpace<I>> pIndexSpace,
                CSRView<const T,const I> lhs) -> std::shared_ptr<LoggedOperator<DefaultSpace<T>>> {
                    CIE_BEGIN_EXCEPTION_TRACING
                        auto pLHSOperator = std::make_shared<CSROperator<I,T,T>>(
                            lhs,
                            pScalarSpace->getThreads());
                        const cie::io::JSONObject preconditionerConfiguration = rConfiguration["preconditioner"];
                        std::optional<std::shared_ptr<LinearOperator<DefaultSpace<T>>>> pMaybePreconditioner;
                        if (!preconditionerConfiguration.is<std::nullptr_t>()) {
                            pMaybePreconditioner = rFactory.make(
                                preconditionerConfiguration,
                                pScalarSpace,
                                pIndexSpace,
                                lhs);
                            CIE_CHECK(
                                pMaybePreconditioner.has_value(),
                                std::format(
                                    "\"{}\" does not name a registered preconditioner",
                                    preconditionerConfiguration["type"].as<std::string>()))
                        }
                        return std::make_shared<ConjugateGradients<DefaultSpace<T>>>(
                            pLHSOperator,
                            pScalarSpace,
                            pMaybePreconditioner.has_value() ? pMaybePreconditioner.value() : nullptr,
                            typename ConjugateGradients<DefaultSpace<T>>::Status {
                                .iterationCount = static_cast<std::size_t>(rConfiguration["max-iterations"].as<double>()),
                                .absoluteResidual = static_cast<T>(rConfiguration["absolute-residual"].as<double>()),
                                .relativeResidual = static_cast<T>(rConfiguration["relative-residual"].as<double>())},
                            static_cast<Verbosity>(rConfiguration["verbosity"].as<int>()));
                    CIE_END_EXCEPTION_TRACING
                };
            rFactory.insert("cg", schema, factory);
        CIE_END_EXCEPTION_TRACING
}


template <class T, class I>
void registerEigenCG(Ref<LinearSolverFactory<
        DefaultSpace<T>,
        DefaultSpace<I>>
    > rFactory) {
        CIE_BEGIN_EXCEPTION_TRACING
            cie::io::JSONObject schema(std::string {R"({"$ref" : "/cie/linalg/cg-eigen"})"});
            const auto factory = [] (
                Ref<const cie::io::JSONObject> rConfiguration,
                std::shared_ptr<DefaultSpace<T>> pScalarSpace,
                std::shared_ptr<DefaultSpace<I>>,
                CSRView<const T,const I> lhs) -> std::shared_ptr<LoggedOperator<DefaultSpace<T>>> {
                    CIE_BEGIN_EXCEPTION_TRACING
                        return std::make_shared<EigenCG<T,I>>(
                            lhs,
                            static_cast<I>(rConfiguration["max-iterations"].as<double>()),
                            static_cast<T>(rConfiguration["relative-residual"].as<double>()),
                            rConfiguration["preconditioner"].as<std::string>(),
                            static_cast<Verbosity>(rConfiguration["verbosity"].as<int>()),
                            pScalarSpace->getThreads());
                    CIE_END_EXCEPTION_TRACING
                };
            rFactory.insert("cg-eigen", schema, factory);
        CIE_END_EXCEPTION_TRACING
}


template <class T, class I>
void registerEigenLLT(Ref<LinearSolverFactory<
        DefaultSpace<T>,
        DefaultSpace<I>>
    > rFactory) {
        CIE_BEGIN_EXCEPTION_TRACING
            cie::io::JSONObject schema(std::string {R"({"$ref" : "/cie/linalg/llt-eigen"})"});
            const auto factory = [] (
                Ref<const cie::io::JSONObject>,
                std::shared_ptr<DefaultSpace<T>> pScalarSpace,
                std::shared_ptr<DefaultSpace<I>>,
                CSRView<const T,const I> lhs) -> std::shared_ptr<LoggedOperator<DefaultSpace<T>>> {
                    CIE_BEGIN_EXCEPTION_TRACING
                        return makeLoggedOperator<EigenLLT<T,I>>(
                            lhs,
                            pScalarSpace->getThreads());
                    CIE_END_EXCEPTION_TRACING
                };
            rFactory.insert("llt-eigen", schema, factory);
        CIE_END_EXCEPTION_TRACING
}


template <class T, class I>
void registerEigenLDLT(Ref<LinearSolverFactory<
        DefaultSpace<T>,
        DefaultSpace<I>>
    > rFactory) {
        CIE_BEGIN_EXCEPTION_TRACING
            cie::io::JSONObject schema(std::string {R"({"$ref" : "/cie/linalg/ldlt-eigen"})"});
            const auto factory = [] (
                Ref<const cie::io::JSONObject>,
                std::shared_ptr<DefaultSpace<T>> pScalarSpace,
                std::shared_ptr<DefaultSpace<I>>,
                CSRView<const T,const I> lhs) -> std::shared_ptr<LoggedOperator<DefaultSpace<T>>> {
                    CIE_BEGIN_EXCEPTION_TRACING
                        return makeLoggedOperator<EigenLDLT<T,I>>(
                            lhs,
                            pScalarSpace->getThreads());
                    CIE_END_EXCEPTION_TRACING
                };
            rFactory.insert("ldlt-eigen", schema, factory);
        CIE_END_EXCEPTION_TRACING
}


#ifdef CIE_ENABLE_SYCL


template <class T, class I>
void registerInverseDiagonal(Ref<LinearSolverFactory<
        SYCLSpace<T>,
        SYCLSpace<I>>
    > rFactory) {
        CIE_BEGIN_EXCEPTION_TRACING
            cie::io::JSONObject schema(std::string {R"({"$ref" : "/cie/linalg/inverse-diagonal-operator-sycl"})"});
            const auto factory = [] (
                Ref<const cie::io::JSONObject>,
                std::shared_ptr<SYCLSpace<T>> pScalarSpace,
                std::shared_ptr<SYCLSpace<I>>,
                CSRView<const T,const I> lhs) -> std::shared_ptr<LoggedOperator<SYCLSpace<T>>> {
                    CIE_BEGIN_EXCEPTION_TRACING
                        return makeLoggedOperator<DiagonalOperator<SYCLSpace<T>>>(
                            makeDiagonalOperator<T,I,T>(
                                lhs,
                                pScalarSpace));
                    CIE_END_EXCEPTION_TRACING
                };
            rFactory.insert("inverse-diagonal-operator-sycl", schema, factory);
        CIE_END_EXCEPTION_TRACING
}


template <class T, class I>
void registerJacobi(Ref<LinearSolverFactory<
        SYCLSpace<T>,
        SYCLSpace<I>>
    > rFactory) {
        CIE_BEGIN_EXCEPTION_TRACING
            cie::io::JSONObject schema(std::string {R"({"$ref" : "/cie/linalg/jacobi-sycl"})"});
            const auto factory = [] (
                Ref<const cie::io::JSONObject> rConfiguration,
                std::shared_ptr<SYCLSpace<T>> pScalarSpace,
                std::shared_ptr<SYCLSpace<I>>,
                CSRView<const T,const I> lhs) -> std::shared_ptr<LoggedOperator<SYCLSpace<T>>> {
                    CIE_BEGIN_EXCEPTION_TRACING
                        auto pLHSOperator = std::make_shared<SYCLCSROperator<I,T,T>>(
                            lhs,
                            pScalarSpace);
                        auto pDiagonalOperator = std::make_shared<DiagonalOperator<SYCLSpace<T>>>(
                            makeDiagonalOperator<T,I,T>(
                                lhs,
                                pScalarSpace));
                        return makeLoggedOperator<JacobiOperator<SYCLSpace<T>>>(
                            pScalarSpace,
                            lhs.rowCount(),
                            pLHSOperator,
                            pDiagonalOperator,
                            rConfiguration["iterations"].as<double>(),
                            rConfiguration["relaxation"].as<double>());
                    CIE_END_EXCEPTION_TRACING
                };
            rFactory.insert("jacobi-sycl", schema, factory);
        CIE_END_EXCEPTION_TRACING
}


template <class T, class I>
void registerCG(Ref<LinearSolverFactory<
        SYCLSpace<T>,
        SYCLSpace<I>>
    > rFactory) {
        CIE_BEGIN_EXCEPTION_TRACING
            cie::io::JSONObject schema(std::string {R"({"$ref" : "/cie/linalg/cg-sycl"})"});
            const auto factory = [&rFactory] (
                Ref<const cie::io::JSONObject> rConfiguration,
                std::shared_ptr<SYCLSpace<T>> pScalarSpace,
                std::shared_ptr<SYCLSpace<I>> pIndexSpace,
                CSRView<const T,const I> lhs) -> std::shared_ptr<LoggedOperator<SYCLSpace<T>>> {
                    CIE_BEGIN_EXCEPTION_TRACING
                        auto pLHSOperator = std::make_shared<SYCLCSROperator<I,T,T>>(
                            lhs,
                            pScalarSpace);
                        const cie::io::JSONObject preconditionerConfiguration = rConfiguration["preconditioner"];
                        std::optional<std::shared_ptr<LinearOperator<SYCLSpace<T>>>> pMaybePreconditioner;
                        if (!preconditionerConfiguration.is<std::nullptr_t>()) {
                            pMaybePreconditioner = rFactory.make(
                                preconditionerConfiguration,
                                pScalarSpace,
                                pIndexSpace,
                                lhs);
                            CIE_CHECK(
                                pMaybePreconditioner.has_value(),
                                std::format(
                                    "\"{}\" does not name a registered preconditioner",
                                    preconditionerConfiguration["type"].as<std::string>()))
                        }
                        return std::make_shared<ConjugateGradients<SYCLSpace<T>>>(
                            pLHSOperator,
                            pScalarSpace,
                            pMaybePreconditioner.has_value() ? pMaybePreconditioner.value() : nullptr,
                            typename ConjugateGradients<SYCLSpace<T>>::Status {
                                .iterationCount = static_cast<std::size_t>(rConfiguration["max-iterations"].as<double>()),
                                .absoluteResidual = static_cast<T>(rConfiguration["absolute-residual"].as<double>()),
                                .relativeResidual = static_cast<T>(rConfiguration["relative-residual"].as<double>())},
                            static_cast<Verbosity>(rConfiguration["verbosity"].as<int>()));
                    CIE_END_EXCEPTION_TRACING
                };
            rFactory.insert("cg-sycl", schema, factory);
        CIE_END_EXCEPTION_TRACING
}


template <class T, class I>
void registerDefaultToSYCLCast(Ref<LinearSolverFactory<
        DefaultSpace<T>,
        DefaultSpace<I>>
    > rFactory,
    std::string_view reference) {
        CIE_BEGIN_EXCEPTION_TRACING
            cie::io::JSONObject schema(
                std::format(
                    "{}\"$ref\" : \"cie/linalg/{}\"{}",
                    '{', reference, '}'));
            const auto factory = [] (
                Ref<const cie::io::JSONObject> rConfiguration,
                std::shared_ptr<DefaultSpace<T>> pScalarSpace,
                std::shared_ptr<DefaultSpace<I>> pIndexSpace,
                CSRView<const T,const I> lhs) -> std::shared_ptr<LoggedOperator<DefaultSpace<T>>> {
                    CIE_BEGIN_EXCEPTION_TRACING
                        // Construct target spaces.
                        using TargetScalarSpace = SYCLSpace<T>;
                        using TargetIndexSpace = SYCLSpace<I>;
                        const auto pQueue = std::make_shared<sycl::queue>(sycl::default_selector_v);;
                        const auto pTargetScalarSpace = std::make_shared<TargetScalarSpace>(pQueue);
                        const auto pTargetIndexSpace = std::make_shared<TargetIndexSpace>(pQueue);

                        // Construct cast operator.
                        using TargetFactorySingleton = LinearSolverFactorySingleton<TargetScalarSpace,TargetIndexSpace>;
                        const auto pOp = std::make_shared<CastSolver<
                            DefaultSpace<T>,DefaultSpace<I>,
                            TargetScalarSpace,TargetIndexSpace>>(
                                lhs,
                                pScalarSpace,
                                pIndexSpace,
                                pTargetScalarSpace,
                                pTargetIndexSpace,
                                rConfiguration,
                                TargetFactorySingleton::get());
                        return pOp;
                    CIE_END_EXCEPTION_TRACING
                };
            rFactory.insert(reference, schema, factory);
        CIE_END_EXCEPTION_TRACING
}


#endif


template <LinalgSpaceLike TSS, LinalgSpaceLike TIS>
void LinearSolverFactory<TSS,TIS>::load() {
    if constexpr (std::is_same_v<TSS,DefaultSpace<typename TSS::Value>> && std::is_same_v<TIS,DefaultSpace<typename TIS::Value>>) {
        CIE_BEGIN_EXCEPTION_TRACING
            // Smoothers and preconditioners.
            registerInverseDiagonal(*this);
            registerJacobi(*this);

            #ifdef CIE_ENABLE_SYCL
                registerDefaultToSYCLCast(*this, "inverse-diagonal-operator-sycl");
                registerDefaultToSYCLCast(*this, "jacobi-sycl");
            #endif

            // Iterative solvers.
            registerCG(*this);

            #ifdef CIE_ENABLE_SYCL
                registerDefaultToSYCLCast(*this, "cg-sycl");
            #endif

            // Solvers wrapped from Eigen.
            registerEigenCG(*this);
            registerEigenLLT(*this);
            registerEigenLDLT(*this);
        CIE_END_EXCEPTION_TRACING
    }
    #ifdef CIE_ENABLE_SYCL
    else if constexpr (std::is_same_v<TSS,SYCLSpace<typename TSS::Value>> && std::is_same_v<TIS,SYCLSpace<typename TIS::Value>>) {
        CIE_BEGIN_EXCEPTION_TRACING
            registerInverseDiagonal(*this);
            registerJacobi(*this);
            registerCG(*this);
        CIE_END_EXCEPTION_TRACING
    }
    #endif
    else static_assert(std::is_same_v<TSS,void> && std::is_same_v<TIS,void>, "unsupported LinalgSpace");
}


} // namespace cie::linalg


namespace cie {


#define CIE_INSTANTIATE_LINEAR_SOLVER_FACTORY(TS, T, I)                 \
    template class FactoryRegistry<                                     \
        linalg::LoggedOperator<TS<T>>,                                  \
        std::shared_ptr<TS<T>>,                                         \
        std::shared_ptr<TS<I>>,                                         \
        linalg::CSRView<                                                \
            const T,                                                    \
            const I>>;                                                  \
    template class linalg::LinearSolverFactory<TS<T>,TS<I>>;            \
    template class linalg::LinearSolverFactory<TS<T>,TS<I>>::Singleton< \
        linalg::LinearSolverFactory<TS<T>,TS<I>>>;

CIE_INSTANTIATE_LINEAR_SOLVER_FACTORY(linalg::DefaultSpace, float, int)
CIE_INSTANTIATE_LINEAR_SOLVER_FACTORY(linalg::DefaultSpace, double, int)

#ifdef CIE_ENABLE_SYCL
    CIE_INSTANTIATE_LINEAR_SOLVER_FACTORY(linalg::SYCLSpace, float, int)
    CIE_INSTANTIATE_LINEAR_SOLVER_FACTORY(linalg::SYCLSpace, double, int)
#endif

#undef CIE_INSTANTIATE_LINEAR_SOLVER_FACTORY


} // namespace cie::linalg
