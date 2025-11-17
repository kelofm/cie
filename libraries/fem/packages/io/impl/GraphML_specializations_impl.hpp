// help the language server
#include "packages/io/inc/GraphML_specializations.hpp"

// --- FEM Includes ---
#include "packages/io/inc/GraphML_specializations.hpp"

// --- Utility Includes ---
#include "packages/compile_time/packages/parameter_pack/inc/Match.hpp"

// --- STL Includes ---
#include <string> // std::string
#include <iomanip> // std::scientific, std::setprecision
#include <optional>



namespace cie::fem::io {


inline void GraphML::Deserializer<std::string>::onElementBegin(Ptr<void>,
                                                               std::string_view elementName,
                                                               std::span<GraphML::AttributePair>)
{
    if (elementName != "txt") {
        CIE_THROW(
            Exception,
            "Expecting a <txt> element but got <" << elementName << ">."
        )
    }
}


inline void GraphML::Deserializer<std::string>::onText(Ptr<void> pThis,
                                                       std::string_view data)
{
    Ref<Deserializer> rThis = *static_cast<Ptr<Deserializer>>(pThis);
    rThis.instance() = data;
}


inline void GraphML::Deserializer<std::string>::onElementEnd(Ptr<void> pThis,
                                                             std::string_view elementName)
{
    if (elementName != "txt") {
        CIE_THROW(
            Exception,
            "Expecting a closing tag for <txt> but got one for <" << elementName << ">."
        )
    }
    Ref<Deserializer> rThis = *static_cast<Ptr<Deserializer>>(pThis);
    rThis.template release<Deserializer>(&rThis, elementName);
}


template <concepts::Container TContainer>
requires std::is_same_v<
    decltype(
        std::declval<GraphML::Serializer<typename TContainer::value_type>>()(
            std::declval<Ref<GraphML::XMLElement>>(),
            std::declval<Ref<const typename TContainer::value_type>>())),
    void
>
void GraphML::Serializer<TContainer>::header(Ref<XMLElement> rElement)
{
    XMLElement defaultElement = rElement.addChild("default");
    TContainer instance;
    this->operator()(defaultElement, instance);
}


template <concepts::Container TContainer>
requires std::is_same_v<
    decltype(
        std::declval<GraphML::Serializer<typename TContainer::value_type>>()(
            std::declval<Ref<GraphML::XMLElement>>(),
            std::declval<Ref<const typename TContainer::value_type>>())),
    void
>
void GraphML::Serializer<TContainer>::operator()(Ref<XMLElement> rElement,
                                                 Ref<const TContainer> rInstance)
{
    GraphML::XMLElement subElement = rElement.addChild("list");
    subElement.addAttribute("size", std::to_string(rInstance.size()));

    GraphML::Serializer<typename TContainer::value_type> subSerializer;
    for (const auto& rItem : rInstance) {
        //XMLElement child = subElement.addChild("li");
        //subSerializer(child, rItem);
        subSerializer(subElement, rItem);
    }
}


template <concepts::Container TContainer>
requires std::is_same_v<
    decltype(
        GraphML::Deserializer<typename TContainer::value_type>::onText(
            std::declval<Ptr<void>>(),
            std::declval<std::string_view>())),
    void
>
void GraphML::Deserializer<TContainer>::onElementBegin(Ptr<void> pThis,
                                                       std::string_view elementName,
                                                       std::span<GraphML::AttributePair> attributes)
{
    CIE_BEGIN_EXCEPTION_TRACING
    using SubDeserializer = GraphML::Deserializer<typename TContainer::value_type>;
    Ref<Deserializer> rThis = *static_cast<Ptr<Deserializer>>(pThis);

    if (!rThis.maybeIt.has_value()) {
        // Parse container size.
        const auto itSize = std::find_if(attributes.begin(),
                                         attributes.end(),
                                         [](const auto pair) {
                                            return pair.first == "size";
                                         });

        if (itSize == attributes.end()) {
            CIE_THROW(
                Exception,
                "Missing attribute: \"size\" on <" << elementName << ">."
            )
        }

        long long size = 0l;
        auto [pEnd, error] = std::from_chars(itSize->second.data(),
                                             itSize->second.data() + itSize->second.size(),
                                             size);
        if (error != std::errc {} || pEnd != itSize->second.data() + itSize->second.size()) {
            CIE_THROW(
                Exception,
                "Failed to convert \"" << itSize->second << "\" to an integer while parsing GraphML."
            )
        }

        if (size < 0l) {
            CIE_THROW(Exception, "Found a negative size " << size << " while parsing GraphML.")
        }

        // Fill the container.
        const std::size_t arraySize = size;

        if (arraySize) {
            if constexpr (concepts::detail::HasResize<TContainer,std::size_t>) {
                rThis.instance().resize(arraySize);
            } else if constexpr (concepts::StaticContainer<TContainer>) {
                if (rThis.instance().size() != arraySize) {
                CIE_THROW(
                    OutOfRangeException,
                    "Expecting a container of size " << rThis.instance().size() << " "
                    << "but parsed a <" << elementName << "> element with \"size\" " << arraySize <<"."
                )
            }
            } else {
                static_assert(std::is_same_v<TContainer,void>, "unsupported container");
            }
        }

        rThis.maybeIt = rThis.instance().begin();
        rThis.maybeItEnd = rThis.instance().end();
    } /*if !rThis.maybeIt.has_value()*/ else if (rThis.maybeIt.value() != rThis.maybeItEnd.value()) {
        Ptr<SubDeserializer> pSubDeserializer = SubDeserializer::make(*rThis.maybeIt.value(), rThis.sax(), elementName);
        rThis.sax().push({
            pSubDeserializer,
            SubDeserializer::onElementBegin,
            SubDeserializer::onText,
            SubDeserializer::onElementEnd
        });
        ++rThis.maybeIt.value();
        SubDeserializer::onElementBegin(pSubDeserializer, elementName, attributes);
    } else {
        std::stringstream message;
        message << "Unexpected child of <list> at index "
                << std::distance(rThis.instance().begin(), rThis.maybeIt.value())
                << ".";
        CIE_THROW(
            OutOfRangeException,
            message.view()
        )
    }
    CIE_END_EXCEPTION_TRACING
}


template <concepts::Container TContainer>
requires std::is_same_v<
    decltype(
        GraphML::Deserializer<typename TContainer::value_type>::onText(
            std::declval<Ptr<void>>(),
            std::declval<std::string_view>())),
    void
>
void GraphML::Deserializer<TContainer>::onText(Ptr<void>,
                                               std::string_view)
{
    CIE_THROW(
        Exception,
        "Found unexpected text data while parsing a container in GraphML."
    )
}


template <concepts::Container TContainer>
requires std::is_same_v<
    decltype(
        GraphML::Deserializer<typename TContainer::value_type>::onText(
            std::declval<Ptr<void>>(),
            std::declval<std::string_view>())),
    void
>
void GraphML::Deserializer<TContainer>::onElementEnd(Ptr<void> pThis,
                                                     std::string_view elementName)
{
    if (elementName != "list") {
        CIE_THROW(
            Exception,
            "Expecting a closing tag for <list> but got one for <" << elementName << ">."
        )
    }

    Ref<Deserializer> rThis = *static_cast<Ptr<Deserializer>>(pThis);

    if (rThis.maybeIt.has_value() && rThis.maybeIt.value() != rThis.maybeItEnd.value()) {
        std::stringstream message;
        message << "Unexpected end of list while parsing GraphML. Expecting " << rThis.instance().size() << " "
                << "items, but parsed " << std::distance(rThis.instance().begin(), rThis.maybeIt.value()) << ".";
        CIE_THROW(Exception, message.view())
    }

    rThis.template release<Deserializer>(&rThis, elementName);
}


template <concepts::Container TContainer>
requires std::is_arithmetic_v<typename TContainer::value_type>
GraphML::Serializer<TContainer>::Serializer() noexcept
{
    this->setFormat(tags::Text::flags());
}


template <concepts::Container TContainer>
requires std::is_arithmetic_v<typename TContainer::value_type>
void GraphML::Serializer<TContainer>::header(Ref<XMLElement> rElement)
{
    XMLElement defaultElement = rElement.addChild("default");
    TContainer instance;
    this->operator()(defaultElement, instance);
}


template <concepts::Container TContainer>
requires std::is_arithmetic_v<typename TContainer::value_type>
void GraphML::Serializer<TContainer>::operator()(Ref<XMLElement> rElement,
                                                 Ref<const TContainer> rInstance)
{
    CIE_BEGIN_EXCEPTION_TRACING
    GraphML::XMLElement subElement = rElement.addChild("arr");
    subElement.addAttribute("size", std::to_string(rInstance.size()));

    if ((this->_format & tags::Text::flags()) != tags::Flags()) {
        std::stringstream data;

        if constexpr (std::is_same_v<typename TContainer::value_type,float>) {
            data << std::scientific << std::setprecision(8);
        } else if constexpr (std::is_same_v<typename TContainer::value_type,double>) {
            data << std::scientific << std::setprecision(16);
        }

        for (auto item : rInstance) {
            data << item << ' ';
        }

        subElement.setValue(data.view());
    } /*if format == text*/ else {
        CIE_THROW(
            NotImplementedException,
            "Writing arithmetc arrays in binary requires an implementation for Base64 encoding."
        )
    }

    CIE_END_EXCEPTION_TRACING
}


template <concepts::Container TContainer>
requires std::is_arithmetic_v<typename TContainer::value_type>
void GraphML::Serializer<TContainer>::setFormat(tags::Flags format)
{
    if (tags::Flags() == (format & (tags::Text::flags() | tags::Binary::flags()))) {
        CIE_THROW(
            Exception,
            "Invalid arithmetic array format " << format << "."
        )
    }
    _format = format;
}


template <concepts::Container TContainer>
requires std::is_arithmetic_v<typename TContainer::value_type>
GraphML::Deserializer<TContainer>::Deserializer(Ref<TContainer> rInstance,
                                                Ref<GraphML::SAXHandler> rSAX)
    : GraphML::DeserializerBase<TContainer>(rInstance, rSAX)
{
    this->setFormat(tags::Text::flags());
}


template <concepts::Container TContainer>
requires std::is_arithmetic_v<typename TContainer::value_type>
void GraphML::Deserializer<TContainer>::onElementBegin(Ptr<void> pThis,
                                                       std::string_view elementName,
                                                       std::span<GraphML::AttributePair> attributes)
{
    CIE_BEGIN_EXCEPTION_TRACING
    Ref<Deserializer> rThis = *static_cast<Ptr<Deserializer>>(pThis);

    // Parse container size.
    const auto itSize = std::find_if(attributes.begin(),
                                     attributes.end(),
                                     [](const auto pair) {
                                        return pair.first == "size";
                                     });

    if (itSize == attributes.end()) {
        CIE_THROW(
            Exception,
            "Missing attribute: \"size\" on <" << elementName << ">."
        )
    }

    long long size = 0l;
    auto [pEnd, error] = std::from_chars(itSize->second.data(),
                                         itSize->second.data() + itSize->second.size(),
                                         size);
    if (error != std::errc {} || pEnd != itSize->second.data() + itSize->second.size()) {
        CIE_THROW(
            Exception,
            "Failed to convert \"" << itSize->second << "\" to an integer while parsing GraphML."
        )
    }

    if (size < 0l) {
        CIE_THROW(Exception, "Found a negative size " << size << " while parsing GraphML.")
    }

    // Fill the container.
    if (size) {
        if constexpr (concepts::detail::HasResize<TContainer,std::size_t>) {
            rThis.instance().resize(size);
        } else if constexpr (concepts::StaticContainer<TContainer>) {
            if (rThis.instance().size() != static_cast<std::size_t>(size)) {
                CIE_THROW(
                    OutOfRangeException,
                    "Expecting a container of size " << rThis.instance().size() << " "
                    << "but parsed a <" << elementName << "> element with \"size\" " << size <<"."
                )
            }
        } else {
            static_assert(std::is_same_v<TContainer,void>, "unsupported container");
        }
    }

    CIE_END_EXCEPTION_TRACING
}


template <concepts::Container TContainer>
requires std::is_arithmetic_v<typename TContainer::value_type>
void GraphML::Deserializer<TContainer>::onText(Ptr<void> pThis,
                                               std::string_view data)
{
    Ref<Deserializer> rThis = *static_cast<Ptr<Deserializer>>(pThis);
    using TValue = typename TContainer::value_type;

    rThis._maybeParsedCount = 0ul;
    Ref<std::size_t> rParsedCount = rThis._maybeParsedCount.value();
    auto itTarget = rThis.instance().begin() + rParsedCount;

    if ((rThis._format & tags::Text::flags()) != tags::Flags()) {
        // Make sure the last character is a delimiter.
        if (!data.empty()) {
            if (data.back() != ' ') {
                CIE_THROW(
                    Exception,
                    "Expecting text block representing a container to end with whitespace ' ', " <<
                    "but the last character is '" << data.back() << ' ';
                )
            }
        } else {
            return;
        }

        // Parse until the last delimiter.
        CIE_BEGIN_EXCEPTION_TRACING
        const char* itSource = data.begin();
        char* itSourceEnd = nullptr;

        while (itSource != itSourceEnd && itTarget != rThis.instance().end()) {
            if constexpr (std::is_same_v<TValue,float>) {
                char* itParsedEnd = nullptr;
                *itTarget = std::strtof(itSource, &itParsedEnd);
                itSource = itParsedEnd;
            } else if constexpr (std::is_same_v<TValue,double>) {
                char* itParsedEnd = nullptr;
                *itTarget = std::strtod(itSource, &itParsedEnd);
                itSource = itParsedEnd;
            } else {
                static_assert(std::is_same_v<TValue,void>, "Unsupported arithmetic type.");
            }

            ++itTarget;
        }
        CIE_END_EXCEPTION_TRACING
    } /*if text*/ else {
        CIE_THROW(
            NotImplementedException,
            "Writing arithmetc arrays in binary requires an implementation for Base64 encoding."
        )
    }

    rParsedCount = std::distance(rThis.instance().begin(), itTarget);
}


template <concepts::Container TContainer>
requires std::is_arithmetic_v<typename TContainer::value_type>
void GraphML::Deserializer<TContainer>::onElementEnd(Ptr<void> pThis,
                                                     std::string_view elementName)
{
    if (elementName != "arr") {
        CIE_THROW(
            Exception,
            "Expecting a closing tag for <arr> but got one for <" << elementName << ">."
        )
    }

    Ref<Deserializer> rThis = *static_cast<Ptr<Deserializer>>(pThis);
    if (!rThis._maybeParsedCount.has_value()) {
        if (!rThis.instance().empty()) {
            CIE_THROW(
                Exception,
                "Expecting a text block for <" << elementName << "> but got none."
            )
        }
    } else if (rThis._maybeParsedCount.value() != rThis.instance().size()) {
        CIE_THROW(
            Exception,
            "Expecting " << rThis.instance().size() << " "
            << " items in the text block of <" << elementName << "> "
            << "but got " << rThis._maybeParsedCount.value() << ".";
        )
    }
    rThis.template release<Deserializer>(&rThis, elementName);
}


template <concepts::Container TContainer>
requires std::is_arithmetic_v<typename TContainer::value_type>
void GraphML::Deserializer<TContainer>::setFormat(tags::Flags format)
{
    if (tags::Flags() == (format & (tags::Text::flags() | tags::Binary::flags()))) {
        CIE_THROW(
            Exception,
            "Invalid arithmetic array format " << format << "."
        )
    }
    _format = format;
}


} // namespace cie::fem::io
