#pragma once

// --- Utility Includes ---
#include "packages/io/inc/json.hpp"

// --- STL Includes ---
#include <concepts>
#include <memory>
#include <functional>
#include <optional>
#include <string_view>
#include <vector>
#include <tuple>
#include <mutex>


namespace cie {


template <class T, class ...TArgs>
class FactoryRegistry {
public:
    using FactoryFunctor = std::function<std::shared_ptr<T>(Ref<const cie::io::JSONObject>,TArgs...)>;

    using Value = T;

    using ArgTuple = std::tuple<TArgs...>;

    FactoryRegistry();

    virtual ~FactoryRegistry();

    [[nodiscard]] std::optional<std::shared_ptr<T>> make(
        Ref<const cie::io::JSONObject> rConfiguration,
        TArgs...) const;

    [[nodiscard]] virtual std::optional<std::shared_ptr<T>> make(
        std::string_view key,
        Ref<const cie::io::JSONObject> rConfiguration,
        TArgs... rArgs) const;

    [[nodiscard]] virtual std::optional<cie::io::JSONObject> makeSchema(std::string_view key) const;

    virtual void insert(
        std::string_view key,
        Ref<const cie::io::JSONObject> rSchema,
        Ref<const FactoryFunctor> rFactory);

    virtual void erase(std::string_view key);

    [[nodiscard]] virtual std::vector<std::string_view> keys() const;

    virtual void load();

    template <class TRegistry = FactoryRegistry<T,TArgs...>>
    class Singleton {
    public:
        [[nodiscard]] static Ref<TRegistry> get();

        static void clear();

    private:
        Singleton() = delete;

        ~Singleton() = delete;

        static inline std::optional<TRegistry> _maybeRegistry = {};

        static inline std::mutex _mutex = {};
    }; // class Singleton

private:
    struct Impl;
    std::unique_ptr<Impl> _pImpl;
}; // class FactoryRegistry


} // namespace cie
