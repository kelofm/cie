// --- External Includes ---
#include "libxml/tree.h"

// --- FEM Includes ---
#include "packages/io/inc/GraphML.hpp"

// --- STL Includes ---
#include <algorithm> // std::copy
#include <cstring> // std::memset
#include <stack> // std::stack
#include <iostream> // std::cin


namespace cie::fem::io {


struct GraphML::XMLElement::Impl
{
    xmlNodePtr _pNode;
}; // struct GraphML::XMLElement::Impl


GraphML::XMLElement::XMLElement()
    : XMLElement(nullptr)
{
}


GraphML::XMLElement::XMLElement(void* pWrapped)
    : _pImpl(new Impl {static_cast<xmlNode*>(pWrapped)})
{
}


GraphML::XMLElement::XMLElement(XMLElement&&) noexcept = default;


GraphML::XMLElement::~XMLElement() = default;


void GraphML::XMLElement::addAttribute(std::string_view key,
                                        std::string_view value)
{
    std::vector<xmlChar> keyCopy, valueCopy;

    keyCopy.resize(key.size() + 1);
    std::copy(key.begin(), key.end(), keyCopy.begin());
    keyCopy.back() = '\0';

    valueCopy.resize(value.size() + 1);
    std::copy(value.begin(), value.end(), valueCopy.begin());
    valueCopy.back() = '\0';

    xmlNewProp(_pImpl->_pNode, keyCopy.data(), valueCopy.data());
}


void GraphML::XMLElement::setValue(std::string_view value)
{
    std::vector<xmlChar> buffer(value.size() + 1);
    std::copy(value.begin(), value.end(), buffer.begin());
    buffer.back() = '\0';
    xmlNodeSetContent(_pImpl->_pNode, buffer.data());
}


GraphML::XMLElement GraphML::XMLElement::addChild(std::string_view name)
{
    XMLElement child;

    std::vector<xmlChar> buffer(name.size() +1);
    std::copy(name.begin(), name.end(), buffer.begin());
    buffer.back() = '\0';

    child._pImpl->_pNode = xmlNewChild(_pImpl->_pNode,
                                       nullptr,
                                       buffer.data(),
                                       nullptr);

    return child;
}


struct GraphML::Input::Impl
{
    Ptr<std::istream> _pStream;
}; // struct GraphML::Input::Impl


GraphML::Input::Input()
    : Input(std::cin)
{
}


GraphML::Input::Input(Input&&) noexcept = default;


GraphML::Input::Input(Ref<std::istream> rStream)
    : _pImpl(new Impl {&rStream})
{
}


GraphML::Input::~Input() = default;


Ref<std::istream> GraphML::Input::stream() noexcept
{
    return *_pImpl->_pStream;
}


struct GraphML::SAXHandler::Impl
{
    std::istream* _pStream;

    xmlSAXHandler _wrapped;

    xmlParserCtxtPtr _pParserContext;

    std::stack<State> _stateStack;

    static XMLStringView makeView(const xmlChar* pBegin) noexcept
    {
        const xmlChar* pEnd = pBegin;
        while (*pEnd != '\0') {++pEnd;}
        return {pBegin, pEnd};
    }

    static void onElementBegin(void* pContext,
                               const xmlChar* pLocalName,
                               [[maybe_unused]] const xmlChar* pPrefix,
                               [[maybe_unused]] const xmlChar* pURI,
                               [[maybe_unused]] int namespaceCount,
                               [[maybe_unused]] const xmlChar** pNamespaces,
                               int attributeCount,
                               [[maybe_unused]] int defaultAttributeCount,
                               const xmlChar** pAttributes)
    {
        Ref<SAXHandler> rSAX = *static_cast<Ptr<SAXHandler>>(pContext);
        if (rSAX._pImpl->_stateStack.empty()) {
            CIE_THROW(
                OutOfRangeException,
                "state stack empty while parsing element '" << pLocalName << "'")
        }

        XMLStringView tagName = makeView(pLocalName);
        std::vector<std::pair<XMLStringView,XMLStringView>> attributes(attributeCount);
        for (int iAttribute=0ul; iAttribute<attributeCount; ++iAttribute) {
            attributes[iAttribute].second = makeView(pAttributes[iAttribute]);
        }

        Ref<State> rState = rSAX._pImpl->_stateStack.top();
        std::get<0>(rState)(std::get<3>(rState),
                            rSAX,
                            tagName,
                            {attributes.data(), attributes.data() + attributeCount});
    }

    static void onText(void* pContext,
                       const xmlChar* pBegin,
                       int size)
    {
        Ref<SAXHandler> rSAX = *static_cast<Ptr<SAXHandler>>(pContext);

        if (rSAX._pImpl->_stateStack.empty()) {
            std::string characters('\0', size);
            std::copy(pBegin, pBegin + size, characters.begin());
            CIE_THROW(
                OutOfRangeException,
                "state stack empty while parsing text '" << characters << "'")
        }

        Ref<State> rState = rSAX._pImpl->_stateStack.top();
        std::get<1>(rState)(std::get<3>(rState),
                            rSAX,
                            {pBegin, pBegin + size});
    }

    static void onElementEnd(void* pContext,
                             const xmlChar* pLocalName,
                             [[maybe_unused]] const xmlChar* pPrefix,
                             [[maybe_unused]] const xmlChar* pURI)
    {
        Ref<SAXHandler> rSAX = *static_cast<Ptr<SAXHandler>>(pContext);

        if (rSAX._pImpl->_stateStack.empty()) {
            CIE_THROW(
                OutOfRangeException,
                "state stack empty while parsing the end of element '" << pLocalName << "'")
        }

        Ref<State> rState = rSAX._pImpl->_stateStack.top();
        std::get<2>(rState)(std::get<3>(rState),
                            rSAX,
                            makeView(pLocalName));
        rSAX._pImpl->_stateStack.pop();
    }
}; // struct GraphML::SAXHandler::Impl


GraphML::SAXHandler::SAXHandler(std::istream& rStream)
    : _pImpl(new Impl)
{
    _pImpl->_pStream = &rStream;
    std::memset(&_pImpl->_wrapped, 0, sizeof(xmlSAXHandler));

    _pImpl->_pParserContext = xmlCreatePushParserCtxt(&_pImpl->_wrapped,
                                                      static_cast<void*>(this),
                                                      nullptr,
                                                      0,
                                                      nullptr);
    if (!_pImpl->_pParserContext) CIE_THROW(Exception, "call to xmlCreatePushParserCtxt failed");

    _pImpl->_wrapped.startElementNs = &Impl::onElementBegin;
    _pImpl->_wrapped.characters = &Impl::onText;
    _pImpl->_wrapped.endElementNs = &Impl::onElementEnd;
}


GraphML::SAXHandler::~SAXHandler()
{
    xmlFreeParserCtxt(_pImpl->_pParserContext);
    //xmlCleanupParser();
}


void GraphML::SAXHandler::push(State state)
{
    _pImpl->_stateStack.push(std::move(state));
}


void GraphML::SAXHandler::parse(std::size_t bufferSize)
{
    CIE_BEGIN_EXCEPTION_TRACING

    std::string buffer;
    buffer.resize(bufferSize);

    Ref<std::istream> rStream = *_pImpl->_pStream;
    while (!rStream.eof()) {
        if (rStream.bad())
            CIE_THROW(Exception, "input stream in invalid state while SAX parsing GraphML")

        const auto readSize = rStream.readsome(buffer.data(), bufferSize);
        const auto error = xmlParseChunk(_pImpl->_pParserContext, buffer.data(), readSize, 0);

        if (error)
            CIE_THROW(Exception, "'xmlParseChunk' failed with error code " << error)
    } // while !rStream.eof()

    // Tell the XML parser that the last chunk was parsed.
    xmlParseChunk(_pImpl->_pParserContext, buffer.data(), 0, 1);

    // Fetch validation flags from the parser.
    if (!_pImpl->_pParserContext->wellFormed)
        CIE_THROW(Exception, "input GraphML is not a well formed XML document")

    CIE_END_EXCEPTION_TRACING
}


struct GraphML::Output::Impl
{
    std::filesystem::path _outputPath;

    xmlDocPtr _pDocument;

    xmlNodePtr _pRoot;
}; // struct GraphML::Output::Impl


GraphML::Output::Output()
    : Output("-")
{
}


GraphML::Output::Output(Ref<const std::filesystem::path> rOutputPath)
    : _pImpl(new Impl {._outputPath = rOutputPath,
                       ._pDocument = nullptr,
                       ._pRoot = nullptr})
{
    constexpr xmlChar xmlStandard[] = "1.0";
    _pImpl->_pDocument = xmlNewDoc(xmlStandard);

    constexpr xmlChar rootName[] = "graphml";
    _pImpl->_pRoot = xmlNewNode(nullptr, rootName);
    xmlDocSetRootElement(_pImpl->_pDocument, _pImpl->_pRoot);
}


GraphML::Output::~Output()
{
    xmlFreeDoc(_pImpl->_pDocument);
};


GraphML::XMLElement GraphML::Output::root()
{
    return XMLElement(_pImpl->_pRoot);
}


void GraphML::Output::write()
{
    const auto outputPath = _pImpl->_outputPath.string();
    std::vector<char> buffer(outputPath.size() + 1);
    std::copy(outputPath.begin(), outputPath.end(), buffer.begin());
    buffer.back() = '\0';

    xmlSaveFile(buffer.data(), _pImpl->_pDocument);
}


} // namespace cie::fem::io
