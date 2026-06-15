// --- Linalg Includes ---
#include "packages/constrained_solvers/inc/ConstrainedSolverFactory.hpp"
#include "packages/constrained_solvers/inc/PenaltyConstrainedOperator.hpp"
#include "packages/constrained_solvers/inc/ALMConstrainedOperator.hpp"
#include "packages/solvers/inc/DefaultSpace.hpp"

// --- Utility Includes ---
#include "packages/macros/inc/exceptions.hpp"
#include "packages/io/inc/json.hpp"
#include "packages/registry/impl/FactoryRegistry_impl.hpp"


namespace cie::linalg {


template <class T, class I>
void registerPenaltySolver(Ref<ConstrainedSolverFactory<
        DefaultSpace<T>,
        DefaultSpace<I>>
    > rFactory) {
        CIE_BEGIN_EXCEPTION_TRACING
            cie::io::JSONObject schema(std::string {R"({"$ref" : "/cie/linalg/penalty"})"});
            const auto factory = [] (
                Ref<const cie::io::JSONObject> rConfiguration,
                std::shared_ptr<DefaultSpace<T>> pScalarSpace,
                std::shared_ptr<DefaultSpace<I>> pIndexSpace,
                CSRView<const T,const I> constraintGradients,
                typename DefaultSpace<T>::ConstVectorView constraintGaps,
                CSRView<const T,const I> unconstrainedLHS,
                OptionalRef<LinearSolverFactory<DefaultSpace<T>,DefaultSpace<I>>> rMaybeLinearSolverFactory)
                    -> std::shared_ptr<ConstrainedOperator<DefaultSpace<T>,DefaultSpace<I>>> {
                        CIE_BEGIN_EXCEPTION_TRACING
                            return std::make_shared<PenaltyConstrainedOperator<T,I>>(
                                rConfiguration,
                                pScalarSpace,
                                pIndexSpace,
                                constraintGradients,
                                constraintGaps,
                                unconstrainedLHS,
                                rMaybeLinearSolverFactory);
                        CIE_END_EXCEPTION_TRACING
                };
            rFactory.insert("penalty", schema, factory);
        CIE_END_EXCEPTION_TRACING
}


template <class T, class I>
void registerAugmentedLagrangeSolver(Ref<ConstrainedSolverFactory<
        DefaultSpace<T>,
        DefaultSpace<I>>
    > rFactory) {
        CIE_BEGIN_EXCEPTION_TRACING
            cie::io::JSONObject schema(std::string {R"({"$ref" : "/cie/linalg/augmented-lagrange"})"});
            const auto factory = [] (
                Ref<const cie::io::JSONObject> rConfiguration,
                std::shared_ptr<DefaultSpace<T>> pScalarSpace,
                std::shared_ptr<DefaultSpace<I>> pIndexSpace,
                CSRView<const T,const I> constraintGradients,
                typename DefaultSpace<T>::ConstVectorView constraintGaps,
                CSRView<const T,const I> unconstrainedLHS,
                OptionalRef<LinearSolverFactory<DefaultSpace<T>,DefaultSpace<I>>> rMaybeLinearSolverFactory)
                    -> std::shared_ptr<ConstrainedOperator<DefaultSpace<T>,DefaultSpace<I>>> {
                        CIE_BEGIN_EXCEPTION_TRACING
                            return std::make_shared<ALMConstrainedOperator<T,I>>(
                                rConfiguration,
                                pScalarSpace,
                                pIndexSpace,
                                constraintGradients,
                                constraintGaps,
                                unconstrainedLHS,
                                rMaybeLinearSolverFactory);
                        CIE_END_EXCEPTION_TRACING
                };
            rFactory.insert("augmented-lagrange", schema, factory);
        CIE_END_EXCEPTION_TRACING
}


template <LinalgSpaceLike TSS, LinalgSpaceLike TIS>
void ConstrainedSolverFactory<TSS,TIS>::load() {
    if constexpr (std::is_same_v<TSS,DefaultSpace<typename TSS::Value>> && std::is_same_v<TIS,DefaultSpace<typename TIS::Value>>) {
        CIE_BEGIN_EXCEPTION_TRACING
            registerPenaltySolver(*this);
            registerAugmentedLagrangeSolver(*this);
        CIE_END_EXCEPTION_TRACING
    }
    else static_assert(std::is_same_v<TSS,void>, "unsupported LinalgSpace");
}


} // namespace cie::linalg


namespace cie {


#define CIE_INSTANTIATE_LINEAR_SOLVER_FACTORY(TS, T, I)                         \
    template class FactoryRegistry<linalg::ConstrainedOperator<TS<T>,TS<I>>,    \
        std::shared_ptr<TS<T>>,                                                 \
        std::shared_ptr<TS<I>>,                                                 \
        linalg::CSRView<const T,const I>,                                       \
        typename TS<T>::ConstVectorView,                                        \
        linalg::CSRView<const T,const I>,                                       \
        OptionalRef<linalg::LinearSolverFactory<TS<T>,TS<I>>>>;                 \
    template class linalg::ConstrainedSolverFactory<TS<T>,TS<I>>;               \
    template class linalg::ConstrainedSolverFactory<TS<T>,TS<I>>::Singleton<    \
        linalg::ConstrainedSolverFactory<TS<T>,TS<I>>>;

CIE_INSTANTIATE_LINEAR_SOLVER_FACTORY(linalg::DefaultSpace, float, int)
CIE_INSTANTIATE_LINEAR_SOLVER_FACTORY(linalg::DefaultSpace, double, int)

#undef CIE_INSTANTIATE_LINEAR_SOLVER_FACTORY


} // namespace cie::linalg
