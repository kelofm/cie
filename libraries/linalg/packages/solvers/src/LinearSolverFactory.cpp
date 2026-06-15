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
                        auto pLHSOperator = std::make_shared<CSROperator<I,T,T>>(
                            lhs,
                            pScalarSpace->getThreads());
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
            cie::io::JSONObject schema(std::string {R"({"$ref" : "/cie/linalg/cg"})"});
            const auto factory = [&rFactory] (
                Ref<const cie::io::JSONObject> rConfiguration,
                std::shared_ptr<SYCLSpace<T>> pScalarSpace,
                std::shared_ptr<SYCLSpace<I>> pIndexSpace,
                CSRView<const T,const I> lhs) -> std::shared_ptr<LoggedOperator<SYCLSpace<T>>> {
                    CIE_BEGIN_EXCEPTION_TRACING
                        auto pLHSOperator = std::make_shared<CSROperator<I,T,T>>(
                            lhs,
                            pScalarSpace->getThreads());
                        cie::io::JSONObject preconditionerConfiguration = rConfiguration["preconditioner"];
                        auto pMaybePreconditioner = rFactory.make(
                            preconditionerConfiguration,
                            pScalarSpace,
                            pIndexSpace,
                            lhs);
                        CIE_CHECK(
                            pMaybePreconditioner.has_value(),
                            std::format(
                                "\"{}\" does not name a registered preconditioner",
                                preconditionerConfiguration["type"].as<std::string>()))
                        return std::make_shared<ConjugateGradients<SYCLSpace<T>>>(
                            pLHSOperator,
                            pScalarSpace,
                            pMaybePreconditioner.value(),
                            typename ConjugateGradients<SYCLSpace<T>>::Statistics {
                                .iterationCount = static_cast<std::size_t>(rConfiguration["max-iterations"].as<double>()),
                                .absoluteResidual = static_cast<T>(rConfiguration["absolute-residual"].as<double>()),
                                .relativeResidual = static_cast<T>(rConfiguration["relative-residual"].as<double>())},
                            static_cast<Verbosity>(rConfiguration["verbosity"].as<int>()));
                    CIE_END_EXCEPTION_TRACING
                };
            rFactory.insert("cg-sycl", schema, factory);
        CIE_END_EXCEPTION_TRACING
}


#endif


template <LinalgSpaceLike TSS, LinalgSpaceLike TIS>
void LinearSolverFactory<TSS,TIS>::load() {
    if constexpr (std::is_same_v<TSS,DefaultSpace<typename TSS::Value>> && std::is_same_v<TIS,DefaultSpace<typename TIS::Value>>) {
        CIE_BEGIN_EXCEPTION_TRACING
            registerInverseDiagonal(*this);
            registerJacobi(*this);
            registerCG(*this);
            registerEigenCG(*this);
            registerEigenLLT(*this);
            registerEigenLDLT(*this);
        CIE_END_EXCEPTION_TRACING
    }
    #ifdef CIE_ENABLE_SYCL
    else if (std::is_same_v<TSS,SYCLSpace<typename TSS::Value>> && std::is_same_v<TIS,SYCLSpace<typename TIS::Value>>) {
        CIE_BEGIN_EXCEPTION_TRACING
            registerInverseDiagonal(*this);
            registerJacobi(*this);
            registerCG(*this);
        CIE_END_EXCEPTION_TRACING
    }
    #endif
    else static_assert(std::is_same_v<TSS,void>, "unsupported LinalgSpace");
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
