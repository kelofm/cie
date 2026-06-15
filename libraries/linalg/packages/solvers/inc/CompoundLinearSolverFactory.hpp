#pragma once

// --- Linalg Includes ---
#include "packages/solvers/inc/LinearSolverFactory.hpp"

// --- Utility Includes ---
#include "packages/compile_time/packages/parameter_pack/inc/Match.hpp"

// --- STL Includes ---
#include <tuple>


namespace cie::linalg {


template <class ...TFactories>
class CompoundLinearSolverFactory {
public:
    CompoundLinearSolverFactory() noexcept = default;

    CompoundLinearSolverFactory(Ref<const TFactories>... rFactories);

    template <class T, class I>
    requires (ct::Match<T>::template Any<typename TFactories::Scalar...> && ct::Match<I>::template Any<typename TFactories::Index...>)
    std::optional<std::variant<
        std::shared_ptr<typename TFactories::Value>...
    >> make(
        Ref<const cie::io::JSONObject> rConfiguration,
        std::tuple<std::shared_ptr<typename TFactories::ScalarSpace>...> pScalarSpace,
        std::tuple<std::shared_ptr<typename TFactories::IndexSpace>...> pIndexSpace,
        CSRView<const T, const I> lhs,
        bool allowAmbiguous = false) const;

    cie::io::JSONObject makeSchema(Ref<const std::string> rName) const;

    std::vector<std::string> keys() const;

private:
    std::tuple<const TFactories*...> _pFactories;
}; // class CompoundLinearSolverFactory


template <class ...TFactories>
CompoundLinearSolverFactory<TFactories...>
makeCompoundLinearSolverFactory(const TFactories& ... rFactories) {
    return CompoundLinearSolverFactory<TFactories...>(rFactories...);
}


} // namespace cie::linalg

#include "packages/solvers/impl/CompoundLinearSolverFactory_impl.hpp"
