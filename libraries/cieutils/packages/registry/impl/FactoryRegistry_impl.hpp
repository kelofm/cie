#pragma once

// --- Utility Includes ---
#include "packages/registry/inc/FactoryRegistry.hpp"
#include "packages/macros/inc/exceptions.hpp"
#include "packages/macros/inc/checks.hpp"

// --- STL Includes ---
#include <unordered_map>
#include <ranges>
#include <format>


namespace cie {


template <class T, class ...TArgs>
struct FactoryRegistry<T,TArgs...>::Impl {
    struct Entry {
        typename FactoryRegistry<T,TArgs...>::FactoryFunctor factory;
        cie::io::JSONObject schema;
        cie::io::JSONSchema validator;
    };

    std::unordered_map<
        std::string,
        Entry
    > map;
}; // struct FactoryRegistry::Impl


template <class T, class ...TArgs>
FactoryRegistry<T,TArgs...>::FactoryRegistry()
    : _pImpl(new Impl)
{}


template <class T, class ...TArgs>
FactoryRegistry<T,TArgs...>::~FactoryRegistry() = default;


template <class T, class ...TArgs>
std::optional<std::shared_ptr<T>> FactoryRegistry<T,TArgs...>::make(
    Ref<const cie::io::JSONObject> rConfiguration,
    TArgs... rArgs) const {
        CIE_BEGIN_EXCEPTION_TRACING
            CIE_CHECK(
                rConfiguration.hasKey("type"),
                "FactoryRegistry expects the provided configuration "
                    << "to specify a solver's name as \"type\", but got "
                    << rConfiguration)
            return this->make(
                rConfiguration["type"].as<std::string>(),
                rConfiguration,
                std::forward<TArgs>(rArgs)...);
        CIE_END_EXCEPTION_TRACING
}


template <class T, class ...TArgs>
std::optional<std::shared_ptr<T>> FactoryRegistry<T,TArgs...>::make(
    std::string_view key,
    Ref<const cie::io::JSONObject> rConfiguration,
    TArgs... rArgs) const {
        const std::string buffer(key.begin(), key.end());
        const auto it = _pImpl->map.find(buffer);
        if (it == _pImpl->map.end()) return {};
        CIE_BEGIN_EXCEPTION_TRACING
            it->second.validator.validate(rConfiguration);
        CIE_END_EXCEPTION_TRACING
        CIE_BEGIN_EXCEPTION_TRACING
            return it->second.factory(
                rConfiguration,
                std::forward<TArgs>(rArgs)...);
        CIE_END_EXCEPTION_TRACING
}


template <class T, class ...TArgs>
std::optional<cie::io::JSONObject> FactoryRegistry<T,TArgs...>::makeSchema(std::string_view key) const {
    const std::string buffer(key.begin(), key.end());
    const auto it = _pImpl->map.find(buffer);
    if (it == _pImpl->map.end()) return {};
    else return it->second.schema;
}


template <class T, class ...TArgs>
void FactoryRegistry<T,TArgs...>::insert(
    std::string_view key,
    Ref<const cie::io::JSONObject> rSchema,
    Ref<const FactoryFunctor> rFactory) {
        CIE_BEGIN_EXCEPTION_TRACING
            typename Impl::Entry entry {
                    .factory = rFactory,
                    .schema = rSchema,
                    .validator = cie::io::JSONSchema(rSchema)};
            const std::string buffer(key.begin(), key.end());
            const bool inserted = _pImpl->map.emplace(
                buffer,
                std::move(entry)).second;
            CIE_CHECK(
                inserted,
                std::format(
                    "\"{}\" is already registered",
                    key))
        CIE_END_EXCEPTION_TRACING
}



template <class T, class ...TArgs>
void FactoryRegistry<T,TArgs...>::erase(std::string_view key) {
    CIE_BEGIN_EXCEPTION_TRACING
        const std::string buffer(key.begin(), key.end());
        const auto it = _pImpl->map.find(buffer);
        if (it != _pImpl->map.end())
            _pImpl->map.erase(it);
    CIE_END_EXCEPTION_TRACING
}


template <class T, class ...TArgs>
std::vector<std::string_view> FactoryRegistry<T,TArgs...>::keys() const {
    CIE_BEGIN_EXCEPTION_TRACING
        std::vector<std::string_view> output;
        const auto keyRange = std::ranges::views::keys(_pImpl->map);
        output.insert(
            output.end(),
            keyRange.begin(),
            keyRange.end());
        return output;
    CIE_END_EXCEPTION_TRACING
}


template <class T, class ...TArgs>
void FactoryRegistry<T,TArgs...>::load() {}


template <class T, class ...TArgs>
template <class TT>
Ref<TT> FactoryRegistry<T,TArgs...>::Singleton<TT>::get() {
    static_assert(std::derived_from<TT,FactoryRegistry<T,TArgs...>>);
    std::scoped_lock<std::mutex> lock(Singleton::_mutex);
    if (!Singleton::_maybeRegistry.has_value()) {
        Singleton::_maybeRegistry.emplace();
        Singleton::_maybeRegistry.value().load();
    }
    return Singleton::_maybeRegistry.value();
}


template <class T, class ...TArgs>
template <class TT>
void FactoryRegistry<T,TArgs...>::Singleton<TT>::clear() {
    std::scoped_lock<std::mutex> lock(Singleton::_mutex);
    _maybeRegistry.reset();
}


} // namespace cie
