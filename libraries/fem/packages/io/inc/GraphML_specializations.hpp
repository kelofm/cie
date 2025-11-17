#pragma once

// --- FEM Includes --
#include "packages/io/inc/GraphML.hpp"

// --- Utiltiy Includes ---
#include "packages/compile_time/packages/concepts/inc/container_concepts.hpp"
#include "packages/types/inc/tags.hpp"

// --- STL Includes ---
#include <optional>


namespace cie::fem::io {


template <>
struct GraphML::Serializer<void> {};


template <>
struct GraphML::Deserializer<void>
    : public GraphML::DeserializerBase<void>
{
    using GraphML::DeserializerBase<void>::DeserializerBase;

    static void onElementBegin(Ptr<void>,
                               std::string_view,
                               std::span<GraphML::AttributePair>) noexcept
    {}

    static void onText(Ptr<void>,
                       std::string_view)
    {
        CIE_THROW(
            Exception,
            "Unexpected text data while parsing a void element in GraphML."
        )
    }

    static void onElementEnd(Ptr<void> pThis,
                             std::string_view)
    {
        delete static_cast<Ptr<Deserializer>>(pThis);
    }
}; // struct Deserializer<void>


template <>
struct GraphML::Serializer<std::string>
{
    void header(Ref<XMLElement> rElement)
    {
        GraphML::XMLElement child = rElement.addChild("default");
        Serializer()(child, "");
    }

    void operator()(Ref<XMLElement> rElement,
                    Ref<const std::string> rInstance)
    {
        GraphML::XMLElement child = rElement.addChild("txt");
        child.setValue(rInstance);
    }
}; // struct GraphML::Serializer<std::string<arithmetic>>


template <>
struct GraphML::Deserializer<std::string>
    : public GraphML::DeserializerBase<std::string>
{
    using GraphML::DeserializerBase<std::string>::DeserializerBase;

    static void onElementBegin(Ptr<void> pThis,
                               std::string_view name,
                               std::span<GraphML::AttributePair> attributes);

    static void onText(Ptr<void> pThis,
                       std::string_view data);

    static void onElementEnd(Ptr<void> pThis,
                             std::string_view name);
}; // struct Deserializer<std::string>


template <concepts::Container TContainer>
requires std::is_arithmetic_v<typename TContainer::value_type>
struct GraphML::Serializer<TContainer>
{
    Serializer() noexcept;

    void header(Ref<XMLElement>);

    void operator()(Ref<XMLElement>,
                    Ref<const TContainer> rInstance);

    void setFormat(tags::Flags format);

private:
    tags::Flags _format;
}; // struct GraphML::Serializer<TContainer<arithmetic>>


template <concepts::Container TContainer>
requires std::is_arithmetic_v<typename TContainer::value_type>
struct GraphML::Deserializer<TContainer>
    : public GraphML::DeserializerBase<TContainer>
{
    Deserializer(Ref<TContainer> rInstance,
                 Ref<GraphML::SAXHandler> rSAX);

    static void onElementBegin(Ptr<void> pThis,
                               std::string_view name,
                               std::span<GraphML::AttributePair> attributes);

    static void onText(Ptr<void> pThis,
                       std::string_view data);

    static void onElementEnd(Ptr<void> pThis,
                             std::string_view name);

    void setFormat(tags::Flags format);

private:
    std::optional<std::size_t> _maybeParsedCount;

    tags::Flags _format;
}; // struct GraphML::Deserializer<TContainer<arithmetic>>


template <concepts::Container TContainer>
requires std::is_same_v<
    decltype(
        std::declval<GraphML::Serializer<typename TContainer::value_type>>()(
            std::declval<Ref<GraphML::XMLElement>>(),
            std::declval<Ref<const typename TContainer::value_type>>())),
    void
>
struct GraphML::Serializer<TContainer>
{
    void header(Ref<XMLElement>);

    void operator()(Ref<XMLElement>,
                    Ref<const TContainer> rInstance);
}; // struct GraphML::Serializer<TContainer<serializable>>


template <concepts::Container TContainer>
requires std::is_same_v<
    decltype(
        GraphML::Deserializer<typename TContainer::value_type>::onText(
            std::declval<Ptr<void>>(),
            std::declval<std::string_view>())),
    void
>
struct GraphML::Deserializer<TContainer>
    : public GraphML::DeserializerBase<TContainer>
{
    using GraphML::DeserializerBase<TContainer>::DeserializerBase;

    static void onElementBegin(Ptr<void> pThis,
                               std::string_view name,
                               std::span<GraphML::AttributePair> attributes);

    static void onText(Ptr<void> pThis,
                       std::string_view data);

    static void onElementEnd(Ptr<void> pThis,
                             std::string_view name);

    std::optional<typename TContainer::iterator> maybeIt, maybeItEnd;
}; // struct GraphML::Deserializer<TContainer<deserializable>>


} // namespace cie::fem::io

#include "packages/io/impl/GraphML_specializations_impl.hpp"
