#pragma once

// help the language server
#include "packages/io/inc/GraphML.hpp"

// --- Utility Includes ---
#include "packages/compile_time/packages/parameter_pack/inc/Match.hpp"

// --- STL Includes ---
#include <charconv> // std::from_chars
#include <iostream>


namespace cie::fem::io {


template <class T>
Ptr<GraphML::Deserializer<T>> GraphML::DeserializerBase<T>::make(typename VoidSafe<T,std::nullptr_t>::Ref rInstance,
                                                                 Ref<SAXHandler> rSAX,
                                                                 [[maybe_unused]] std::string_view elementName)
{
    return new GraphML::Deserializer<T>(rInstance, rSAX);
}


template <class T>
template <class TDerived>
void GraphML::DeserializerBase<T>::release(Ptr<TDerived> pThis,
                                           [[maybe_unused]] std::string_view elementName) noexcept
{
    //std::cout << "release " << "\033[0;31m" << elementName << "\033[0;37m"
    //          << " " << typeid(T).name()
    //          << std::endl;
    delete pThis;
}


template <class T>
GraphML::DeserializerBase<T>::DeserializerBase(typename VoidSafe<T,std::nullptr_t>::Ref rInstance,
                                               Ref<SAXHandler> rSAX) noexcept
    : _pInstance(&rInstance),
      _pSAX(&rSAX)
{
}


template <class T>
typename VoidSafe<T,std::nullptr_t>::Ref GraphML::DeserializerBase<T>::instance() noexcept
{
    if constexpr (std::is_same_v<T,void>) {
        return _pInstance;
    } else {
        return *_pInstance;
    }
}


template <class T>
Ref<GraphML::SAXHandler> GraphML::DeserializerBase<T>::sax() noexcept
{
    return *_pSAX;
}


namespace detail {


struct GraphMLNoOpElement {};


struct GraphMLKeyElement {};


template <class TGraph>
struct GraphMLElement {
    Ptr<TGraph> pGraph;
};


template <class TGraph>
struct GraphMLRoot {
    GraphMLElement<TGraph> graphML;
};


} // namespace detail


template <class TGraph>
struct GraphML::Deserializer<detail::GraphMLElement<TGraph>>
    : public GraphML::DeserializerBase<detail::GraphMLElement<TGraph>>
{
    using GraphML::DeserializerBase<detail::GraphMLElement<TGraph>>::DeserializerBase;

    static void onElementBegin(Ptr<void> pThis,
                               std::string_view elementName,
                               std::span<GraphML::AttributePair>)
    {
        Ref<Deserializer> rThis = *static_cast<Ptr<Deserializer>>(pThis);

        if (elementName == "graph") {
            using SubDeserializer = GraphML::Deserializer<TGraph>;
            rThis.sax().push({
                SubDeserializer::make(*rThis.instance().pGraph, rThis.sax(), elementName),
                SubDeserializer::onElementBegin,
                SubDeserializer::onText,
                SubDeserializer::onElementEnd
            });
        } else {
            CIE_THROW(
                Exception,
                "Expecting a <graph> element, but got <" << elementName << ">."
            )
        }
    }

    static void onText(Ptr<void>,
                       std::string_view)
    {
        CIE_THROW(
            Exception,
            "Unexpected text data on a <graphml> element."
        )
    }

    static void onElementEnd(Ptr<void> pThis,
                             [[maybe_unused]] std::string_view elementName)
    {
        if (elementName != "graphml") {
            CIE_THROW(
                Exception,
                "Expecting a </graphml> end but got </" << elementName << ">."
            )
        }
        Ref<Deserializer> rThis = *static_cast<Ptr<Deserializer>>(pThis);
        rThis.template release<Deserializer>(&rThis, elementName);
    }
}; // struct GraphMLDeserializer


template <class TGraph>
struct GraphML::Deserializer<detail::GraphMLRoot<TGraph>>
    : public GraphML::DeserializerBase<detail::GraphMLRoot<TGraph>>
{
    using GraphML::DeserializerBase<detail::GraphMLRoot<TGraph>>::DeserializerBase;

    using GraphML::DeserializerBase<detail::GraphMLRoot<TGraph>>::release;

    static void onElementBegin(Ptr<void> pThis,
                               std::string_view elementName,
                               std::span<GraphML::AttributePair> attributes)
    {
        Ref<Deserializer> rThis = *static_cast<Ptr<Deserializer>>(pThis);

        if (elementName == "graphml") {
            if (!attributes.empty()) {
                CIE_THROW(
                    Exception,
                    "Expecting no attributes for <graphml>, but got " << attributes.size() << " of them."
                )
            }

            using SubDeserializer = GraphML::Deserializer<detail::GraphMLElement<TGraph>>;
            rThis.sax().push({
                SubDeserializer::make(rThis.instance().graphML, rThis.sax(), elementName),
                SubDeserializer::onElementBegin,
                SubDeserializer::onText,
                SubDeserializer::onElementEnd
            });
        } else {
            CIE_THROW(
                Exception,
                "Expecting a <graphml> element, but got <" << elementName << ">."
            )
        }
    }

    static void onText(Ptr<void>,
                       std::string_view)
    {
        CIE_THROW(
            Exception,
            "Unexpected text data outside a <graphml> element."
        )
    }

    static void onElementEnd(Ptr<void> pThis,
                             std::string_view elementName) noexcept
    {
        Ref<Deserializer> rThis = *static_cast<Ptr<Deserializer>>(pThis);
        rThis.template release<Deserializer>(&rThis, elementName);
    }
}; // struct Deserializer<detail::GraphMLRoot>


template <>
struct GraphML::Deserializer<detail::GraphMLNoOpElement>
    : public GraphML::DeserializerBase<detail::GraphMLNoOpElement>
{
    using GraphML::DeserializerBase<detail::GraphMLNoOpElement>::DeserializerBase;

    static void onElementBegin(Ptr<void> pThis,
                               std::string_view elementName,
                               std::span<GraphML::AttributePair>) noexcept
    {
        Ref<Deserializer> rThis = *static_cast<Ptr<Deserializer>>(pThis);
        using SubDeserializer = Deserializer;
        rThis.sax().push({
            SubDeserializer::make(rThis.instance(), rThis.sax(), elementName),
            SubDeserializer::onElementBegin,
            SubDeserializer::onText,
            SubDeserializer::onElementEnd
        });
    }

    static void onText(Ptr<void>,
                       std::string_view)
    {}

    static void onElementEnd(Ptr<void> pThis,
                             [[maybe_unused]] std::string_view elementName) noexcept
    {
        Ref<Deserializer> rThis = *static_cast<Ptr<Deserializer>>(pThis);
        rThis.template release<Deserializer>(&rThis, elementName);
    }
}; // struct NoOpGraphMLDeserializer


template <>
struct GraphML::Deserializer<detail::GraphMLKeyElement>
    : public GraphML::DeserializerBase<detail::GraphMLKeyElement>
{
    using GraphML::DeserializerBase<detail::GraphMLKeyElement>::DeserializerBase;

    static void onElementBegin(Ptr<void> pThis,
                               std::string_view elementName,
                               std::span<GraphML::AttributePair>)
    {
        Ref<Deserializer> rThis = *static_cast<Ptr<Deserializer>>(pThis);

        if (elementName == "default" || elementName == "desc") {
            using SubDeserializer = GraphML::Deserializer<detail::GraphMLNoOpElement>;
            detail::GraphMLNoOpElement dummy;
            rThis.sax().push({
                SubDeserializer::make(dummy, rThis.sax(), elementName),
                SubDeserializer::onElementBegin,
                SubDeserializer::onText,
                SubDeserializer::onElementEnd
            });
        } else {
            CIE_THROW(
                Exception,
                "Unexpected element <" << elementName << "> while parsing <key>."
            )
        }
    }

    static void onText(Ptr<void>,
                       std::string_view)
    {
        CIE_THROW(
            Exception,
            "Unexpected text data on a <key> element."
        )
    }

    static void onElementEnd(Ptr<void> pThis,
                             [[maybe_unused]] std::string_view elementName)
    {
        if (elementName != "key") {
            CIE_THROW(
                Exception,
                "Expecting a closing tag for <key> but got one for <" << elementName << ">."
            )
        }

        Ref<Deserializer> rThis = *static_cast<Ptr<Deserializer>>(pThis);
        rThis.template release<Deserializer>(&rThis, elementName);
    }
}; // struct GraphMLKeyDeserializer



template <class TVertexData, class TEdgeData, class TGraphData>
void GraphML::Input::operator()(Ref<Graph<TVertexData,TEdgeData,TGraphData>> rGraph)
{
    CIE_BEGIN_EXCEPTION_TRACING

    SAXHandler sax(this->stream());

    using SubDeserializer = GraphML::Deserializer<detail::GraphMLRoot<Graph<TVertexData,TEdgeData,TGraphData>>>;
    detail::GraphMLRoot<Graph<TVertexData,TEdgeData,TGraphData>> root;
    root.graphML.pGraph = &rGraph;
    Ptr<SubDeserializer> pSubDeserializer = SubDeserializer::make(root, sax, "");

    sax.push({
        pSubDeserializer,
        SubDeserializer::onElementBegin,
        SubDeserializer::onText,
        SubDeserializer::onElementEnd
    });

    sax.parse();
    pSubDeserializer->template release<SubDeserializer>(pSubDeserializer, "");

    CIE_END_EXCEPTION_TRACING
}


template <class TVertexData, class TEdgeData, class TGraphData>
void GraphML::Output::operator()(Ref<const Graph<TVertexData,TEdgeData,TGraphData>> rGraph)
{
    CIE_BEGIN_EXCEPTION_TRACING
        XMLElement rootElement = this->root();
        XMLElement graphElement = rootElement.addChild("graph");

        Serializer<Graph<TVertexData,TEdgeData,TGraphData>> serializer;
        serializer.header(graphElement);
        serializer(graphElement, rGraph);
    CIE_END_EXCEPTION_TRACING

    CIE_BEGIN_EXCEPTION_TRACING
        this->write();
    CIE_END_EXCEPTION_TRACING
}


template <class TVertexData, class TEdgeData, class TGraphData>
void GraphML::Serializer<Graph<TVertexData,TEdgeData,TGraphData>>::header(Ref<GraphML::XMLElement> rElement) const
{
    CIE_BEGIN_EXCEPTION_TRACING

    int dataCount = 0;    // <== counts how many non-void properties need to be written
    std::string graphDataID, vertexDataID, edgeDataID;

    // Write graph attributes.
    rElement.addAttribute("edgedefault", "directed");

    // Write graph data header.
    if constexpr (ct::Match<TGraphData>::template None<void,std::monostate>) {
        GraphML::Serializer<TGraphData> serializer;
        graphDataID = std::to_string(dataCount++);

        XMLElement headerElement = rElement.addChild("key");
        headerElement.addAttribute("id", graphDataID);
        headerElement.addAttribute("for", "graph");

        serializer.header(headerElement);
    }

    // Write vertex data header.
    if constexpr (ct::Match<TVertexData>::template None<void,std::monostate>) {
        GraphML::Serializer<TVertexData> serializer;
        vertexDataID = std::to_string(dataCount++);

        XMLElement headerElement = rElement.addChild("key");
        headerElement.addAttribute("id", vertexDataID);
        headerElement.addAttribute("for", "node");

        serializer.header(headerElement);
    }

    // Write edge data header.
    if constexpr (ct::Match<TEdgeData>::template None<void,std::monostate>) {
        GraphML::Serializer<TEdgeData> serializer;
        edgeDataID = std::to_string(dataCount++);

        XMLElement headerElement = rElement.addChild("key");
        headerElement.addAttribute("id", edgeDataID);
        headerElement.addAttribute("for", "edge");
        serializer.header(headerElement);
    }

    CIE_END_EXCEPTION_TRACING
}


template <class TVertexData, class TEdgeData, class TGraphData>
void GraphML::Serializer<Graph<TVertexData,TEdgeData,TGraphData>>::operator()(Ref<GraphML::XMLElement> rElement,
                                                                              Ref<const Graph<TVertexData,TEdgeData,TGraphData>> rGraph) const
{
    CIE_BEGIN_EXCEPTION_TRACING

    // Find keys.
    int dataCount = 0;    // <== counts how many non-void properties need to be written
    std::string graphDataID, vertexDataID, edgeDataID;

    if constexpr (ct::Match<TGraphData>::template None<void,std::monostate>) {
        graphDataID = std::to_string(dataCount++);
    }
    if constexpr (ct::Match<TVertexData>::template None<void,std::monostate>) {
        vertexDataID = std::to_string(dataCount++);
    }
    if constexpr (ct::Match<TEdgeData>::template None<void,std::monostate>) {
        edgeDataID = std::to_string(dataCount++);
    }

    // Construct serializers.
    GraphML::Serializer<TGraphData> graphDataSerializer;
    GraphML::Serializer<TVertexData> vertexDataSerializer;
    GraphML::Serializer<TEdgeData> edgeDataSerializer;

    // Write graph data.
    if constexpr (ct::Match<TGraphData>::template None<void,std::monostate>) {
        XMLElement dataElement = rElement.addChild("data");
        dataElement.addAttribute("key", graphDataID);
        graphDataSerializer(dataElement, rGraph.data());
    }

    // Write vertices.
    for (const auto& rItem : rGraph.vertices()) {
        XMLElement element = rElement.addChild("node");
        element.addAttribute("id", std::to_string(rItem.id()));

        if constexpr (ct::Match<TVertexData>::template None<void,std::monostate>) {
            XMLElement dataElement = element.addChild("data");
            vertexDataSerializer(dataElement, rItem.data());
            dataElement.addAttribute("key", vertexDataID);
        }
    } // for rItem in rGraph.vertices()

    // Write edges.
    for (const auto& rItem : rGraph.edges()) {
        XMLElement element = rElement.addChild("edge");
        element.addAttribute("id", std::to_string(rItem.id()));
        element.addAttribute("source", std::to_string(rItem.source()));
        element.addAttribute("target", std::to_string(rItem.target()));

        if constexpr (ct::Match<TEdgeData>::template None<void,std::monostate>) {
            XMLElement dataElement = element.addChild("data");
            edgeDataSerializer(dataElement, rItem.data());
            dataElement.addAttribute("key", edgeDataID);
        }
    } // for rItem in rGraph.edges()

    CIE_END_EXCEPTION_TRACING
}


template <class TVertexData, class TEdgeData, class TGraphData>
void GraphML::Deserializer<Graph<TVertexData,TEdgeData,TGraphData>>::onElementBegin(Ptr<void> pThis,
                                                                                    std::string_view elementName,
                                                                                    std::span<GraphML::AttributePair> attributes)
{
    using Value = Graph<TVertexData,TEdgeData,TGraphData>;
    Ref<Deserializer> rThis = *static_cast<Ptr<Deserializer>>(pThis);

    if (elementName == "node") {
        // Parse the ID.
        const auto itID = std::find_if(attributes.begin(),
                                       attributes.end(),
                                       [](auto pair) {return pair.first == "id";});
        if (itID == attributes.end()) {
            CIE_THROW(Exception, "Found a node without an ID while parsing GraphML.")
        }

        long id = 0l;
        auto [pEnd, error] = std::from_chars(itID->second.data(),
                                             itID->second.data() + itID->second.size(),
                                             id);
        if (error != std::errc {} || pEnd != itID->second.data() + itID->second.size()) {
            CIE_THROW(
                Exception,
                "Failed to convert \"" << itID->second << "\" to a node ID while parsing GraphML."
            )
        }

        if (id < 0l) {
            CIE_THROW(Exception, "Found a node with an invalid ID " << id << " while parsing GraphML.")
        }

        // Construct an instance that will be filled in by its deserializer later.
        Ref<typename Value::Vertex> rItem = rThis.instance().insert(typename Value::Vertex(VertexID(id)));

        // Push a deserializer that fills in the constructed object, and push it
        // onto the top of the parsing stack.
        using SubDeserializer = GraphML::Deserializer<typename Value::Vertex>;
        rThis.sax().push({
            SubDeserializer::make(rItem, rThis.sax(), elementName),
            SubDeserializer::onElementBegin,
            SubDeserializer::onText,
            SubDeserializer::onElementEnd
        });
    } /*if elementName == "node"*/ else if (elementName == "edge") {
        long edgeID     = 0l,
             sourceID   = 0l,
             targetID   = 0l;

        // Parse the edge's ID.
        {
            const auto itID = std::find_if(attributes.begin(),
                                           attributes.end(),
                                           [](auto pair) {return pair.first == "id";});
            if (itID == attributes.end()) {
                CIE_THROW(Exception, "Found an edge without an ID while parsing GraphML.")
            }

            auto [pEnd, error] = std::from_chars(itID->second.data(),
                                                 itID->second.data() + itID->second.size(),
                                                 edgeID);
            if (error != std::errc {} || pEnd != itID->second.data() + itID->second.size()) {
                CIE_THROW(
                    Exception,
                    "Failed to convert \"" << itID->second << "\" to an edge ID while parsing GraphML."
                )
            }
        }

        if (edgeID < 0l) {
            CIE_THROW(Exception, "Found an edge with an invalid ID " << edgeID << " while parsing GraphML.")
        }

        // Parse the source vertex' ID.
        {
            const auto itID = std::find_if(attributes.begin(),
                                           attributes.end(),
                                           [](auto pair) {return pair.first == "source";});
            if (itID == attributes.end()) {
                CIE_THROW(
                    Exception,
                    "The source vertex of edge " << edgeID << " is missing its ID."
                )
            }

            auto [pEnd, error] = std::from_chars(itID->second.data(),
                                                 itID->second.data() + itID->second.size(),
                                                 sourceID);
            if (error != std::errc {} || pEnd != itID->second.data() + itID->second.size()) {
                CIE_THROW(
                    Exception,
                       "Failed to convert \"" << itID->second << "\" "
                    << "to the source vertex ID of edge " << edgeID << "."
                )
            }
        }

        if (sourceID < 0l) {
            CIE_THROW(
                Exception,
                   "Found invalid source vertex ID " << sourceID << " "
                << "of edge " << edgeID << " "
                << "while parsing GraphML."
            )
        }

        // Parse the target vertex' ID.
        {
            const auto itID = std::find_if(attributes.begin(),
                                           attributes.end(),
                                           [](auto pair) {return pair.first == "target";});
            if (itID == attributes.end()) {
                CIE_THROW(
                    Exception,
                    "The target vertex of edge " << edgeID << " is missing its ID."
                )
            }

            auto [pEnd, error] = std::from_chars(itID->second.data(),
                                                 itID->second.data() + itID->second.size(),
                                                 targetID);
            if (error != std::errc {} || pEnd != itID->second.data() + itID->second.size()) {
                CIE_THROW(
                    Exception,
                       "Failed to convert \"" << itID->second << "\" "
                    << "to the target vertex ID of edge " << edgeID << "."
                )
            }
        }

        if (targetID < 0l) {
            CIE_THROW(
                Exception,
                   "Found invalid target vertex ID " << targetID << " "
                << "of edge " << edgeID << " "
                << "while parsing GraphML."
            )
        }

        // Construct an instance that will be filled in by its deserializer later.
        typename Value::Edge& rItem = rThis.instance().insert(typename Value::Edge(
            EdgeID(edgeID),
            {
                VertexID(sourceID),
                VertexID(targetID)
            }
        ));

        // Push a deserializer that fills in the constructed object, and push it
        // onto the top of the parsing stack.
        using SubDeserializer = GraphML::Deserializer<typename Value::Edge>;
        rThis.sax().push({
            SubDeserializer::make(rItem, rThis.sax(), elementName),
            SubDeserializer::onElementBegin,
            SubDeserializer::onText,
            SubDeserializer::onElementEnd
        });
    } /*if elementName == "edge"*/ else if (elementName == "key") {
        using SubDeserializer = GraphML::Deserializer<detail::GraphMLKeyElement>;
        detail::GraphMLKeyElement dummy;
        rThis.sax().push({
            SubDeserializer::make(dummy, rThis.sax(), elementName),
            SubDeserializer::onElementBegin,
            SubDeserializer::onText,
            SubDeserializer::onElementEnd
        });
    } /*if elementName == "key"*/ else if (elementName == "data") {
        using SubDeserializer = GraphML::Deserializer<TGraphData>;
        Ptr<SubDeserializer> pSubDeserializer = SubDeserializer::make(rThis.instance().data(), rThis.sax(), elementName);
        rThis.sax().push({
            pSubDeserializer,
            SubDeserializer::onElementBegin,
            SubDeserializer::onText,
            SubDeserializer::onElementEnd
        });
        //SubDeserializer::onElementBegin(pSubDeserializer, elementName, attributes);
    } /*if elementName == "data"*/ else {
        CIE_THROW(
            Exception,
            "Unexpected <" << elementName << "> "
            << "while parsing <graph>."
        )
    }
}


template <class TVertexData, class TEdgeData, class TGraphData>
void GraphML::Deserializer<Graph<TVertexData,TEdgeData,TGraphData>>::onText(Ptr<void>,
                                                                            std::string_view)
{
    CIE_THROW(
        Exception,
        "Found unexpected text data while parsing GraphML."
    )
}


template <class TVertexData, class TEdgeData, class TGraphData>
void GraphML::Deserializer<Graph<TVertexData,TEdgeData,TGraphData>>::onElementEnd(Ptr<void> pThis,
                                                                                  std::string_view elementName)
{
    Ref<Deserializer> rThis = *static_cast<Ptr<Deserializer>>(pThis);

    if (elementName == "node" || elementName == "edge" || elementName == "data") {
        rThis.sax().push({
            pThis,
            Deserializer::onElementBegin,
            Deserializer::onText,
            Deserializer::onElementEnd
        });
    } else if (elementName == "graph") {
        rThis.template release<Deserializer>(&rThis, elementName);
    } else {
        CIE_THROW(
            Exception,
            "Expecting a closing tag for <graph> but got one for <" << elementName << ">."
        )
    }
}


template <class T>
requires std::is_same_v<typename T::ID,VertexID>
void GraphML::Deserializer<T>::onElementBegin(Ptr<void> pThis,
                                              std::string_view elementName,
                                              [[maybe_unused]] std::span<GraphML::AttributePair> attributes)
{
    Ref<Deserializer> rThis = *static_cast<Ptr<Deserializer>>(pThis);

    if (elementName == "data") {
        using SubDeserializer = GraphML::Deserializer<typename T::Data>;
        Ptr<SubDeserializer> pSubDeserializer;
        if constexpr (std::is_same_v<typename T::Data,void>) {
            pSubDeserializer = SubDeserializer::make(nullptr, rThis.sax(), elementName);
        } else {
            pSubDeserializer = SubDeserializer::make(rThis.instance().data(), rThis.sax(), elementName);
        }

        rThis.sax().push({
            pSubDeserializer,
            SubDeserializer::onElementBegin,
            SubDeserializer::onText,
            SubDeserializer::onElementEnd
        });

        //CIE_BEGIN_EXCEPTION_TRACING
        //SubDeserializer::onElementBegin(pSubDeserializer, elementName, attributes);
        //CIE_END_EXCEPTION_TRACING
    } else {
        CIE_THROW(
            Exception,
            "Expecting a \"data\" graphml element on a \"node\", but got \"" << elementName << "\"."
        )
    }
}


template <class T>
requires std::is_same_v<typename T::ID,VertexID>
void GraphML::Deserializer<T>::onText(Ptr<void>,
                                      std::string_view)
{
    CIE_THROW(
        Exception,
        "Found unexpected text data while parsing a node in GraphML."
    )
}


template <class T>
requires std::is_same_v<typename T::ID,VertexID>
void GraphML::Deserializer<T>::onElementEnd(Ptr<void> pThis,
                                            std::string_view elementName) noexcept
{
    Ref<Deserializer> rThis = *static_cast<Ptr<Deserializer>>(pThis);
    rThis.template release<Deserializer>(&rThis, elementName);
}


template <class T>
requires std::is_same_v<typename T::ID,EdgeID>
void GraphML::Deserializer<T>::onElementBegin(Ptr<void> pThis,
                                              std::string_view elementName,
                                              [[maybe_unused]] std::span<GraphML::AttributePair> attributes)
{
    Ref<Deserializer> rThis = *static_cast<Ptr<Deserializer>>(pThis);

    if (elementName == "data") {
        using SubDeserializer = GraphML::Deserializer<typename T::Data>;
        Ptr<SubDeserializer> pSubDeserializer;
        if constexpr (std::is_same_v<typename T::Data,void>) {
            pSubDeserializer = SubDeserializer::make(nullptr, rThis.sax(), elementName);
        } else {
            pSubDeserializer = SubDeserializer::make(rThis.instance().data(), rThis.sax(), elementName);
        }

        rThis.sax().push({
            pSubDeserializer,
            SubDeserializer::onElementBegin,
            SubDeserializer::onText,
            SubDeserializer::onElementEnd
        });

        //CIE_BEGIN_EXCEPTION_TRACING
        //SubDeserializer::onElementBegin(pSubDeserializer, elementName, attributes);
        //CIE_END_EXCEPTION_TRACING
    } else {
        CIE_THROW(
            Exception,
            "Expecting a \"data\" graphml element on an \"edge\", but got \"" << elementName << "\"."
        )
    }
}


template <class T>
requires std::is_same_v<typename T::ID,EdgeID>
void GraphML::Deserializer<T>::onText(Ptr<void>,
                                      std::string_view)
{
    CIE_THROW(
        Exception,
        "Found unexpected text data while parsing an edge in GraphML."
    )
}


template <class T>
requires std::is_same_v<typename T::ID,EdgeID>
void GraphML::Deserializer<T>::onElementEnd(Ptr<void> pThis,
                                            std::string_view elementName) noexcept
{
    Ref<Deserializer> rThis = *static_cast<Ptr<Deserializer>>(pThis);
    rThis.template release<Deserializer>(&rThis, elementName);
}


} // namespace cie::fem::io
