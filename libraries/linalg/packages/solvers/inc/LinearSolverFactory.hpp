#pragma once

// --- Linalg Includes ---
#include "packages/solvers/inc/LinalgSpace.hpp"
#include "packages/solvers/inc/LinearOperator.hpp"
#include "packages/utilities/inc/CSRView.hpp"

// --- Utility Includes ---
#include "packages/io/inc/json.hpp"

// --- STL Includes ---
#include <string>
#include <functional>
#include <unordered_map>
#include <memory>
#include <ranges>
#include <optional>


namespace cie::linalg {


template <LinalgSpaceLike TScalarSpace, LinalgSpaceLike TIndexSpace>
class LinearSolverFactory {
public:
    using ScalarSpace = TScalarSpace;

    using IndexSpace  = TIndexSpace;

    using Scalar = typename ScalarSpace::Value;

    using Index = typename IndexSpace::Value;

    using Solver = LinearOperator<TScalarSpace>;

    using FactoryFunctor = std::function<
        std::shared_ptr<Solver>(
            Ref<cie::io::JSONObject>,
            std::shared_ptr<ScalarSpace>,
            std::shared_ptr<IndexSpace>,
            CSRView<const Scalar,const Index>)>;

    std::optional<std::shared_ptr<Solver>> make(
        Ref<cie::io::JSONObject> rConfiguration,
        std::shared_ptr<ScalarSpace> pScalarSpace,
        std::shared_ptr<IndexSpace> pIndexSpace,
        CSRView<const Scalar,const Index> lhs) const;

    cie::io::JSONObject makeSchema(Ref<const std::string> rName) const;

    void insert(
        Ref<const std::string> rName,
        Ref<const cie::io::JSONObject> rSchema,
        Ref<const FactoryFunctor> rFactory);

    void erase(Ref<const std::string> rName);

    auto keys() const noexcept {
        return std::ranges::views::keys(_registry);
    }

    class Singleton {
    public:
        static Ref<LinearSolverFactory<TScalarSpace,TIndexSpace>> get();

        static void clear();

    private:
        Singleton() = delete;

        ~Singleton() = delete;

        static std::optional<LinearSolverFactory<TScalarSpace,TIndexSpace>> _maybeFactory;
    }; // class Singleton

private:
    struct Entry {
        cie::io::JSONObject schema;
        FactoryFunctor factory;
    }; // struct Entry

    std::unordered_map<
        std::string,
        Entry
    > _registry;
}; // class LinearSolverFactory


} // namespace cie::linalg
