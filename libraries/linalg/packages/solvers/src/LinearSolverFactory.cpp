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
#include "packages/solvers/inc/SYCLCSROperator.hpp"

// --- Utility Includes ---
#include "packages/macros/inc/exceptions.hpp"
#include "packages/io/inc/json.hpp"


namespace cie::linalg {


template <LinalgSpaceLike TSS, LinalgSpaceLike TIS>
std::optional<std::shared_ptr<LinearOperator<TSS>>> LinearSolverFactory<TSS,TIS>::make(
    Ref<cie::io::JSONObject> rConfiguration,
    std::shared_ptr<ScalarSpace> pScalarSpace,
    std::shared_ptr<IndexSpace> pIndexSpace,
    CSRView<const Scalar,const Index> lhs) const {
        CIE_CHECK(
            rConfiguration.hasKey("type") && rConfiguration["type"].is<std::string>(),
            std::format(
                "expecting linear solver configuration to specify a solver's name as \"type\", but got {}",
                [&rConfiguration](){std::stringstream m; m<<rConfiguration;return m.str();}()))

        CIE_BEGIN_EXCEPTION_TRACING
            const std::string& name = rConfiguration["type"].as<std::string>();
            const auto it = _registry.find(name);
            if (it == _registry.end()) return {};
            CIE_BEGIN_EXCEPTION_TRACING
                cie::io::JSONSchema validator(it->second.schema);
                validator.validateAndFillDefaults(rConfiguration);
            CIE_END_EXCEPTION_TRACING
            CIE_BEGIN_EXCEPTION_TRACING
                return it->second.factory(
                    rConfiguration,
                    pScalarSpace,
                    pIndexSpace,
                    lhs);
            CIE_END_EXCEPTION_TRACING
        CIE_END_EXCEPTION_TRACING
}


template <LinalgSpaceLike TSS, LinalgSpaceLike TIS>
cie::io::JSONObject LinearSolverFactory<TSS,TIS>::makeSchema(Ref<const std::string> rName) const {
    const auto it = _registry.find(rName);
    CIE_CHECK(
        it != _registry.end(),
        std::format(
            "\"{}\" does not name a registered linear solver",
            rName))
    return it->second.schema;
}


template <LinalgSpaceLike TSS, LinalgSpaceLike TIS>
void LinearSolverFactory<TSS,TIS>::insert(
    Ref<const std::string> rName,
    Ref<const cie::io::JSONObject> rSchema,
    Ref<const FactoryFunctor> rFactory) {
        CIE_BEGIN_EXCEPTION_TRACING
            const bool emplaced = _registry.emplace(
                rName,
                Entry {
                    .schema = rSchema,
                    .factory = rFactory}).second;
            CIE_CHECK(
                emplaced,
                std::format(
                    "attempt to re-register linear solver \"{}\"",
                    rName))
        CIE_END_EXCEPTION_TRACING
}


template <LinalgSpaceLike TSS, LinalgSpaceLike TIS>
void LinearSolverFactory<TSS,TIS>::erase(Ref<const std::string> rName) {
    const auto it = _registry.find(rName);
    if (it != _registry.end()) _registry.erase(it);
}


template <class T, class I>
void registerInverseDiagonal(Ref<LinearSolverFactory<
        DefaultSpace<T>,
        DefaultSpace<I>>
    > rFactory) {
        CIE_BEGIN_EXCEPTION_TRACING
            cie::io::JSONObject schema(std::string {R"({"$ref" : "/cie/linalg/inverse-diagonal-operator"})"});
            const auto factory = [] (
                Ref<cie::io::JSONObject>,
                std::shared_ptr<DefaultSpace<T>> pScalarSpace,
                std::shared_ptr<DefaultSpace<I>>,
                CSRView<const T,const I> lhs) -> std::shared_ptr<LinearOperator<DefaultSpace<T>>> {
                    CIE_BEGIN_EXCEPTION_TRACING
                        return std::make_shared<DiagonalOperator<DefaultSpace<T>>>(
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
                Ref<cie::io::JSONObject> rConfiguration,
                std::shared_ptr<DefaultSpace<T>> pScalarSpace,
                std::shared_ptr<DefaultSpace<I>>,
                CSRView<const T,const I> lhs) -> std::shared_ptr<LinearOperator<DefaultSpace<T>>> {
                    CIE_BEGIN_EXCEPTION_TRACING
                        auto pLHSOperator = std::make_shared<CSROperator<I,T,T>>(
                            lhs,
                            pScalarSpace->getThreads());
                        auto pDiagonalOperator = std::make_shared<DiagonalOperator<DefaultSpace<T>>>(
                            makeDiagonalOperator<T,I,T>(
                                lhs,
                                pScalarSpace));
                        return std::make_shared<JacobiOperator<DefaultSpace<T>>>(
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
                Ref<cie::io::JSONObject> rConfiguration,
                std::shared_ptr<DefaultSpace<T>> pScalarSpace,
                std::shared_ptr<DefaultSpace<I>> pIndexSpace,
                CSRView<const T,const I> lhs) -> std::shared_ptr<LinearOperator<DefaultSpace<T>>> {
                    CIE_BEGIN_EXCEPTION_TRACING
                        auto pLHSOperator = std::make_shared<CSROperator<I,T,T>>(
                            lhs,
                            pScalarSpace->getThreads());
                        cie::io::JSONObject preconditionerConfiguration = rConfiguration["preconditioner"];
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
                            typename ConjugateGradients<DefaultSpace<T>>::Statistics {
                                .iterationCount = static_cast<std::size_t>(rConfiguration["max-iterations"].as<double>()),
                                .absoluteResidual = static_cast<T>(rConfiguration["absolute-residual"].as<double>()),
                                .relativeResidual = static_cast<T>(rConfiguration["relative-residual"].as<double>())},
                            rConfiguration["verbosity"].as<int>());
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
                Ref<cie::io::JSONObject> rConfiguration,
                std::shared_ptr<DefaultSpace<T>> pScalarSpace,
                std::shared_ptr<DefaultSpace<I>>,
                CSRView<const T,const I> lhs) -> std::shared_ptr<LinearOperator<DefaultSpace<T>>> {
                    CIE_BEGIN_EXCEPTION_TRACING
                        return std::make_shared<EigenCG<T,I>>(
                            lhs,
                            static_cast<I>(rConfiguration["max-iterations"].as<double>()),
                            static_cast<T>(rConfiguration["absolute-residual"].as<double>()),
                            rConfiguration["preconditioner"].as<std::string>(),
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
                Ref<cie::io::JSONObject>,
                std::shared_ptr<DefaultSpace<T>> pScalarSpace,
                std::shared_ptr<DefaultSpace<I>>,
                CSRView<const T,const I> lhs) -> std::shared_ptr<LinearOperator<DefaultSpace<T>>> {
                    CIE_BEGIN_EXCEPTION_TRACING
                        return std::make_shared<EigenLLT<T,I>>(
                            lhs,
                            pScalarSpace->getThreads());
                    CIE_END_EXCEPTION_TRACING
                };
            rFactory.insert("llt-eigen", schema, factory);
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
                Ref<cie::io::JSONObject>,
                std::shared_ptr<SYCLSpace<T>> pScalarSpace,
                std::shared_ptr<SYCLSpace<I>>,
                CSRView<const T,const I> lhs) -> std::shared_ptr<LinearOperator<SYCLSpace<T>>> {
                    CIE_BEGIN_EXCEPTION_TRACING
                        return std::make_shared<DiagonalOperator<SYCLSpace<T>>>(
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
                Ref<cie::io::JSONObject> rConfiguration,
                std::shared_ptr<SYCLSpace<T>> pScalarSpace,
                std::shared_ptr<SYCLSpace<I>>,
                CSRView<const T,const I> lhs) -> std::shared_ptr<LinearOperator<SYCLSpace<T>>> {
                    CIE_BEGIN_EXCEPTION_TRACING
                        auto pLHSOperator = std::make_shared<CSROperator<I,T,T>>(
                            lhs,
                            pScalarSpace->getThreads());
                        auto pDiagonalOperator = std::make_shared<DiagonalOperator<SYCLSpace<T>>>(
                            makeDiagonalOperator<T,I,T>(
                                lhs,
                                pScalarSpace));
                        return std::make_shared<JacobiOperator<SYCLSpace<T>>>(
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
                Ref<cie::io::JSONObject> rConfiguration,
                std::shared_ptr<SYCLSpace<T>> pScalarSpace,
                std::shared_ptr<SYCLSpace<I>> pIndexSpace,
                CSRView<const T,const I> lhs) -> std::shared_ptr<LinearOperator<SYCLSpace<T>>> {
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
                            rConfiguration["verbosity"].as<int>());
                    CIE_END_EXCEPTION_TRACING
                };
            rFactory.insert("cg-sycl", schema, factory);
        CIE_END_EXCEPTION_TRACING
}


#endif


template <LinalgSpaceLike TSS, LinalgSpaceLike TIS>
Ref<LinearSolverFactory<TSS,TIS>> LinearSolverFactory<TSS,TIS>::Singleton::get() {
    if (!_maybeFactory.has_value()) {
        _maybeFactory.emplace();
        Ref<LinearSolverFactory<TSS,TIS>> rFactory = _maybeFactory.value();

        if constexpr (std::is_same_v<TSS,DefaultSpace<typename TSS::Value>> && std::is_same_v<TIS,DefaultSpace<typename TIS::Value>>) {
            CIE_BEGIN_EXCEPTION_TRACING
                registerInverseDiagonal(rFactory);
                registerJacobi(rFactory);
                registerCG(rFactory);
                registerEigenCG(rFactory);
                registerEigenLLT(rFactory);
            CIE_END_EXCEPTION_TRACING
        }
        #ifdef CIE_ENABLE_SYCL
        else if (std::is_same_v<TSS,SYCLSpace<typename TSS::Value>> && std::is_same_v<TIS,SYCLSpace<typename TIS::Value>>) {
            CIE_BEGIN_EXCEPTION_TRACING
                registerInverseDiagonal(rFactory);
                registerJacobi(rFactory);
                registerCG(rFactory);
            CIE_END_EXCEPTION_TRACING
        }
        #endif
        else static_assert(std::is_same_v<TSS,void>, "unsupported LinalgSpace");
    } // if !_maybeFactory.has_value()
    return _maybeFactory.value();
}


template <LinalgSpaceLike TSS, LinalgSpaceLike TIS>
void LinearSolverFactory<TSS,TIS>::Singleton::clear() {
    _maybeFactory.reset();
}


template <LinalgSpaceLike TSS, LinalgSpaceLike TIS>
std::optional<LinearSolverFactory<TSS,TIS>> LinearSolverFactory<TSS,TIS>::Singleton::_maybeFactory = {};


#define CIE_INSTANTIATE_LINEAR_SOLVER_FACTORY(TSpace, TValue, TIndex)   \
    template class LinearSolverFactory<TSpace<TValue>,TSpace<TIndex>>;

CIE_INSTANTIATE_LINEAR_SOLVER_FACTORY(DefaultSpace, float, int)
CIE_INSTANTIATE_LINEAR_SOLVER_FACTORY(DefaultSpace, double, int)

#ifdef CIE_ENABLE_SYCL
    CIE_INSTANTIATE_LINEAR_SOLVER_FACTORY(SYCLSpace, float, int)
    CIE_INSTANTIATE_LINEAR_SOLVER_FACTORY(SYCLSpace, double, int)
#endif

#undef CIE_INSTANTIATE_LINEAR_SOLVER_FACTORY


} // namespace cie::linalg
