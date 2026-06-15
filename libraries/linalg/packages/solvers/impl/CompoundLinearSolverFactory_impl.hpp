#pragma once

// --- Linalg Includes ---
#include "packages/io/inc/json.hpp"
#include "packages/macros/inc/exceptions.hpp"
#include "packages/solvers/inc/CompoundLinearSolverFactory.hpp"
#include "packages/solvers/inc/DefaultSpace.hpp"


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
                            _lhsRowExtents,
                            _lhsColumnIndices,
                            _lhsEntries);
                        auto pMaybeOperator = rFactory.make(
                            rConfiguration,
                            pTargetScalarSpace,
                            pTargetIndexSpace,
                            targetLhs);
                        _pOperator = pMaybeOperator.value();
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

    std::shared_ptr<LinearOperator<TTargetScalarSpace>> _pOperator;

    typename TTargetIndexSpace::Vector _lhsRowExtents;

    typename TTargetIndexSpace::Vector _lhsColumnIndices;

    typename TTargetScalarSpace::Vector _lhsEntries;

    std::shared_ptr<TSourceScalarSpace> _pSourceScalarSpace;

    std::shared_ptr<TTargetScalarSpace> _pTargetScalarSpace;
}; // class CastSolver


template <class ...TFs>
CompoundLinearSolverFactory<TFs...>::CompoundLinearSolverFactory(Ref<const TFs>... rFactories)
    : _pFactories(&rFactories...)
{}


namespace detail {
template <class T, class I, class ...TFs, std::size_t ...Is>
std::optional<std::variant<
    std::shared_ptr<typename TFs::Value>...
>> compoundFactoryExpansion(
    Ref<const cie::io::JSONObject> rConfiguration,
    std::tuple<const TFs* ...> pFactories,
    std::tuple<std::shared_ptr<typename TFs::ScalarSpace> ...> pScalarSpaces,
    std::tuple<std::shared_ptr<typename TFs::IndexSpace> ...> pIndexSpaces,
    CSRView<const T, const I> lhs,
    bool allowAmbiguous,
    std::index_sequence<Is...>) {
        std::optional<std::variant<
            std::shared_ptr<typename TFs::Value>...
        >> pMaybeOutput;
        CIE_BEGIN_EXCEPTION_TRACING
            CIE_CHECK(
                rConfiguration.hasKey("type"),
                "expecting linear solver configuration to specify a solver's name as \"type\", but got " << rConfiguration)
            const std::string name = rConfiguration["type"].as<std::string>();
            const auto visitor = [&] (
                const auto& rFactory,
                auto pScalarSpace,
                auto pIndexSpace) -> void {
                    // Early exit if ambiguity is allowed and a solver
                    // has already been constructed.
                    if (pMaybeOutput.has_value() && allowAmbiguous)
                        return;

                    // Implicitly assume that the LHS matrix comes from DefaultSpace.
                    using TFactory = std::remove_cvref_t<decltype(rFactory)>;
                    using SourceScalarSpace = DefaultSpace<T>;
                    using SourceIndexSpace = DefaultSpace<I>;
                    using TargetScalarSpace = typename TFactory::ScalarSpace;
                    using TargetIndexSpace = typename TFactory::IndexSpace;

                    // Check whether the solver's name is registered in this factory.
                    const auto& rNames = rFactory.keys();
                    if (std::find(rNames.begin(), rNames.end(), name) != rNames.end()) {
                        CIE_CHECK(
                            !pMaybeOutput.has_value(),
                            '"' << name << '"' << " is ambiguous")

                        if constexpr (std::is_same_v<typename TFactory::ScalarSpace,DefaultSpace<T>> && std::is_same_v<typename TFactory::IndexSpace,DefaultSpace<I>>) {
                            pMaybeOutput = rFactory.make(
                                rConfiguration,
                                pScalarSpace,
                                pIndexSpace,
                                lhs);
                        } else {
                            using Op = CastSolver<SourceScalarSpace,SourceIndexSpace,TargetScalarSpace,TargetIndexSpace>;
                            pMaybeOutput = std::make_shared<Op>(
                                lhs,
                                std::make_shared<SourceScalarSpace>(),
                                std::make_shared<SourceIndexSpace>(),
                                pScalarSpace,
                                pIndexSpace,
                                rConfiguration,
                                rFactory);
                        }
                    }
            }; // visitor
            (visitor(
                *std::get<Is>(pFactories),
                std::get<Is>(pScalarSpaces),
                std::get<Is>(pIndexSpaces))
            ,...);
        CIE_END_EXCEPTION_TRACING
        return pMaybeOutput;
}

template <class ...TFs, std::size_t ...Is>
void compoundFactoryKeyExpansion(
    Ref<std::vector<std::string>> rOutput,
    std::tuple<const TFs*...> pFactories,
    std::index_sequence<Is...>) {
        const auto visitor = [&rOutput] (const auto* pFactory) {
            const auto& rNames = pFactory->keys();
            rOutput.insert(
                rOutput.end(),
                rNames.begin(),
                rNames.end());};
        (visitor(std::get<Is>(pFactories)),...);
}
} // namespace detail


template <class ...TFs>
template <class T, class I>
requires (ct::Match<T>::template Any<typename TFs::Scalar...> && ct::Match<I>::template Any<typename TFs::Index...>)
std::optional<std::variant<
    std::shared_ptr<typename TFs::Value>...
>> CompoundLinearSolverFactory<TFs...>::make(
    Ref<const cie::io::JSONObject> rConfiguration,
    std::tuple<std::shared_ptr<typename TFs::ScalarSpace>...> pScalarSpace,
    std::tuple<std::shared_ptr<typename TFs::IndexSpace>...> pIndexSpace,
    CSRView<const T, const I> lhs,
    bool allowAmbiguous) const {
        CIE_BEGIN_EXCEPTION_TRACING
            return detail::compoundFactoryExpansion<T,I,TFs...>(
                rConfiguration,
                _pFactories,
                pScalarSpace,
                pIndexSpace,
                lhs,
                allowAmbiguous,
                std::make_index_sequence<ct::PackSize<TFs...>>());
        CIE_END_EXCEPTION_TRACING
}


template <class ...TFs>
std::vector<std::string> CompoundLinearSolverFactory<TFs...>::keys() const {
    std::vector<std::string> output;
    CIE_BEGIN_EXCEPTION_TRACING
        detail::compoundFactoryKeyExpansion<TFs...>(
            output,
            _pFactories,
            std::make_index_sequence<ct::PackSize<TFs...>>());
    CIE_END_EXCEPTION_TRACING
    return output;
}


} // namespace cie::linalg
