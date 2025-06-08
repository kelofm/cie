// --- External Includes ---
#include "libxml/tree.h"

// --- FEM Includes ---
#include "packages/io/inc/GraphML.hpp"

// --- STL Includes ---
#include <algorithm> // std::copy


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
