#ifndef CIE_FEM_GRAPHML_IMPL_HPP
#define CIE_FEM_GRAPHML_IMPL_HPP

// help the language server
#include "packages/io/inc/GraphML.hpp"

// --- Utility Includes ---
#include "packages/compile_time/packages/parameter_pack/inc/Match.hpp" // Match::None

// --- STL Includes ---
#include <string> // std::string
#include <charconv> // std::from_chars


namespace cie::fem::io {


namespace detail {


template <class TGraph>
struct GraphMLDeserializer
{
    static void onElementBegin(void* pGraph,
                               Ref<GraphML::SAXHandler> rSAX,
                               std::string_view name,
                               std::span<GraphML::AttributePair> attributes)
    {
        if (name == "graphml") {
            if (!attributes.empty()) {
                CIE_THROW(
                    Exception,
                    "Expecting no attributes for \"graphml\", but got " << attributes.size() << " of them."
                )
            }

            rSAX.push({
                GraphMLDeserializer::onElementBegin,
                GraphMLDeserializer::onText,
                GraphMLDeserializer::onElementEnd,
                pGraph
            });

        } else if (name == "graph") {
            using SubDeserializer = GraphML::Deserializer<TGraph>;
            rSAX.push({
                SubDeserializer::onElementBegin,
                SubDeserializer::onText,
                SubDeserializer::onElementEnd,
                pGraph
            });
        } else {
            CIE_THROW(
                Exception,
                "Expecting a \"graphml\" or \"graph\" element, but got \"" << name << "\"."
            )
        }
    }

    static void onText(void*,
                       Ref<GraphML::SAXHandler>,
                       std::string_view)
    {
        CIE_THROW(
            Exception,
            "Unexpected text data on a \"graphml\" element."
        )
    }

    static void onElementEnd(void*,
                             Ref<GraphML::SAXHandler>,
                             std::string_view) noexcept
    {}
}; // struct GraphMLDeserializer


struct NoOpGraphMLDeserializer
{
    static void onElementBegin(void*,
                               Ref<GraphML::SAXHandler>,
                               std::string_view,
                               std::span<GraphML::AttributePair>) noexcept
    {}

    static void onText(void*,
                       Ref<GraphML::SAXHandler>,
                       std::string_view)
    {}

    static void onElementEnd(void*,
                             Ref<GraphML::SAXHandler>,
                             std::string_view) noexcept
    {}
}; // struct NoOpGraphMLDeserializer


struct GraphMLKeyDeserializer
{
    static void onElementBegin(void* pGraph,
                               Ref<GraphML::SAXHandler> rSAX,
                               std::string_view name,
                               std::span<GraphML::AttributePair>) noexcept
    {
        if (name == "default") {
            rSAX.push({
                NoOpGraphMLDeserializer::onElementBegin,
                NoOpGraphMLDeserializer::onText,
                NoOpGraphMLDeserializer::onElementEnd,
                pGraph
            });
        } else {
            CIE_THROW(
                Exception,
                "Unexpected element \"" << name << "\" while parsing \"key\"."
            )
        }
    }

    static void onText(void*,
                       Ref<GraphML::SAXHandler>,
                       std::string_view)
    {
        CIE_THROW(
            Exception,
            "Unexpected text data on a \"key\" element."
        )
    }

    static void onElementEnd(void*,
                             Ref<GraphML::SAXHandler>,
                             std::string_view) noexcept
    {}
}; // struct GraphMLKeyDeserializer


} // namespace detail


template <class TVertexData, class TEdgeData, class TGraphData>
void GraphML::Input::operator()(Ref<Graph<TVertexData,TEdgeData,TGraphData>> rGraph)
{
    CIE_BEGIN_EXCEPTION_TRACING

    SAXHandler sax(this->stream());

    using SubDeserializer = detail::GraphMLDeserializer<Graph<TVertexData,TEdgeData,TGraphData>>;
    sax.push({
        SubDeserializer::onElementBegin,
        SubDeserializer::onText,
        SubDeserializer::onElementEnd,
        &rGraph
    });

    sax.parse();

    CIE_END_EXCEPTION_TRACING
}


template <>
struct GraphML::Serializer<void> {};


template <class TVertexData, class TEdgeData, class TGraphData>
void GraphML::Output::operator()(Ref<const Graph<TVertexData,TEdgeData,TGraphData>> rGraph)
{
    CIE_BEGIN_EXCEPTION_TRACING

    XMLElement rootElement = this->root();
    XMLElement graphElement = rootElement.addChild("graph");

    // Write graph attributes.
    graphElement.addAttribute("edgedefault", "directed");

    int dataCount = 0;    // <== counts how many non-void properties need to be written
    std::string graphDataID, vertexDataID, edgeDataID;

    // Write graph data header.
    GraphML::Serializer<TGraphData> graphDataSerializer;

    if constexpr (ct::Match<TGraphData>::template None<void,std::monostate>) {
        graphDataID = std::to_string(dataCount++);

        XMLElement headerElement = graphElement.addChild("key");
        headerElement.addAttribute("id", graphDataID);
        headerElement.addAttribute("for", "graph");

        graphDataSerializer.header(headerElement);
    }

    // Write vertex data header.
    GraphML::Serializer<TVertexData> vertexDataSerializer;

    if constexpr (ct::Match<TVertexData>::template None<void,std::monostate>) {
        vertexDataID = std::to_string(dataCount++);

        XMLElement headerElement = graphElement.addChild("key");
        headerElement.addAttribute("id", vertexDataID);
        headerElement.addAttribute("for", "node");

        vertexDataSerializer.header(headerElement);
    }

    // Write edge data header.
    GraphML::Serializer<TEdgeData> edgeDataSerializer;

    if constexpr (ct::Match<TEdgeData>::template None<void,std::monostate>) {
        edgeDataID = std::to_string(dataCount++);

        XMLElement headerElement = graphElement.addChild("key");
        headerElement.addAttribute("id", edgeDataID);
        headerElement.addAttribute("for", "edge");
        edgeDataSerializer.header(headerElement);
    }

    // Write graph data.
    if constexpr (ct::Match<TGraphData>::template None<void,std::monostate>) {
        XMLElement dataElement = graphElement.addChild("data");
        dataElement.addAttribute("key", graphDataID);
        graphDataSerializer(dataElement, rGraph.data());
    }

    // Write vertices.
    for (const auto& rItem : rGraph.vertices()) {
        XMLElement element = graphElement.addChild("node");
        element.addAttribute("id", std::to_string(rItem.id()));

        if constexpr (ct::Match<TVertexData>::template None<void,std::monostate>) {
            XMLElement dataElement = element.addChild("data");
            vertexDataSerializer(dataElement, rItem.data());
            dataElement.addAttribute("key", vertexDataID);
        }
    } // for rItem in rGraph.vertices()

    // Write edges.
    for (const auto& rItem : rGraph.edges()) {
        XMLElement element = graphElement.addChild("edge");
        element.addAttribute("id", std::to_string(rItem.id()));
        element.addAttribute("source", std::to_string(rItem.source()));
        element.addAttribute("target", std::to_string(rItem.target()));

        if constexpr (ct::Match<TEdgeData>::template None<void,std::monostate>) {
            XMLElement dataElement = element.addChild("data");
            edgeDataSerializer(dataElement, rItem.data());
            dataElement.addAttribute("key", edgeDataID);
        }
    } // for rItem in rGraph.edges()

    // Dump.
    this->write();

    CIE_END_EXCEPTION_TRACING
}


template <class TVertexData, class TEdgeData, class TGraphData>
void GraphML::Deserializer<Graph<TVertexData,TEdgeData,TGraphData>>::onElementBegin(void* pGraph,
                                                                                    Ref<GraphML::SAXHandler> rSAX,
                                                                                    std::string_view name,
                                                                                    std::span<GraphML::AttributePair> attributes)
{
    Value& rGraph = *static_cast<Value*>(pGraph);

    if (name == "node") {
        // Parse the ID.
        const auto itID = std::find_if(attributes.begin(),
                                       attributes.end(),
                                       [](auto pair) {return pair.first == "id";});
        if (itID == attributes.end()) {
            for (auto [l, r] : attributes) std::cout << l << " : " << r << std::endl;
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
        typename Value::Vertex& rItem = rGraph.insert(typename Value::Vertex(VertexID(id)));

        // Push a deserializer that fills in the constructed object, and push it
        // onto the top of the parsing stack.
        using SubDeserializer = GraphML::Deserializer<typename Value::Vertex>;
        rSAX.push({
            SubDeserializer::onElementBegin,
            SubDeserializer::onText,
            SubDeserializer::onElementEnd,
            &rItem
        });
    } /*if name == "node"*/ else if (name == "edge") {
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
        typename Value::Edge& rItem = rGraph.insert(typename Value::Edge(
            EdgeID(edgeID),
            {
                VertexID(sourceID),
                VertexID(targetID)
            }
        ));

        // Push a deserializer that fills in the constructed object, and push it
        // onto the top of the parsing stack.
        using SubDeserializer = GraphML::Deserializer<typename Value::Edge>;
        rSAX.push({
            SubDeserializer::onElementBegin,
            SubDeserializer::onText,
            SubDeserializer::onElementEnd,
            &rItem
        });
    } /*if name == "edge"*/ else if (name == "key") {
        rSAX.push({
            detail::GraphMLKeyDeserializer::onElementBegin,
            detail::GraphMLKeyDeserializer::onText,
            detail::GraphMLKeyDeserializer::onElementEnd,
            pGraph
        });
    } /*if name == "key"*/ else if (name == "data") {
        using SubDeserializer = GraphML::Deserializer<TGraphData>;
        rSAX.push({
            SubDeserializer::onElementBegin,
            SubDeserializer::onText,
            SubDeserializer::onElementEnd,
            &rGraph.data()
        });
    } /*if name == "data"*/ else if (name == "graph") {

    } /*if name == "graph"*/ else {
        CIE_THROW(
            Exception,
            "Found unknown element type \"" << name << "\" "
            << "while parsing GraphML."
        )
    }
}


template <class TVertexData, class TEdgeData, class TGraphData>
void GraphML::Deserializer<Graph<TVertexData,TEdgeData,TGraphData>>::onText(void*,
                                                                            Ref<GraphML::SAXHandler>,
                                                                            std::string_view)
{
    CIE_THROW(
        Exception,
        "Found unexpected text data while parsing GraphML."
    )
}


template <class TVertexData, class TEdgeData, class TGraphData>
void GraphML::Deserializer<Graph<TVertexData,TEdgeData,TGraphData>>::onElementEnd(void*,
                                                                                  Ref<GraphML::SAXHandler>,
                                                                                  std::string_view) noexcept
{
}


template <class T>
requires std::is_same_v<typename T::ID,VertexID>
void GraphML::Deserializer<T>::onElementBegin(void* pInstance,
                                              Ref<GraphML::SAXHandler> rSAX,
                                              std::string_view name,
                                              [[maybe_unused]] std::span<GraphML::AttributePair> attributes)
{
    if (name == "data") {
        using SubDeserializer = GraphML::Deserializer<typename T::Data>;
        auto pSubInstance = &static_cast<T*>(pInstance)->data();
        rSAX.push({
            SubDeserializer::onElementBegin,
            SubDeserializer::onText,
            SubDeserializer::onElementEnd,
            pSubInstance
        });

        CIE_BEGIN_EXCEPTION_TRACING
        SubDeserializer::onElementBegin(pSubInstance, rSAX, name, attributes);
        CIE_END_EXCEPTION_TRACING
    } else {
        CIE_THROW(
            Exception,
            "Expecting a \"data\" graphml element on a \"node\", but got \"" << name << "\"."
        )
    }
}


template <class T>
requires std::is_same_v<typename T::ID,VertexID>
void GraphML::Deserializer<T>::onText([[maybe_unused]] void*,
                                      [[maybe_unused]] Ref<GraphML::SAXHandler>,
                                      [[maybe_unused]] std::string_view)
{
    CIE_THROW(
        Exception,
        "Found unexpected text data while parsing a node in GraphML."
    )
}


template <class T>
requires std::is_same_v<typename T::ID,VertexID>
void GraphML::Deserializer<T>::onElementEnd([[maybe_unused]] void* pInstance,
                                            [[maybe_unused]] Ref<GraphML::SAXHandler> rSAX,
                                            [[maybe_unused]] std::string_view name) noexcept
{}


template <class T>
requires std::is_same_v<typename T::ID,EdgeID>
void GraphML::Deserializer<T>::onElementBegin([[maybe_unused]] void* pInstance,
                                              [[maybe_unused]] Ref<GraphML::SAXHandler> rSAX,
                                              [[maybe_unused]] std::string_view name,
                                              [[maybe_unused]] std::span<GraphML::AttributePair> attributes)
{
    if (name == "data") {
        using SubDeserializer = GraphML::Deserializer<typename T::Data>;
        auto pSubInstance = &static_cast<T*>(pInstance)->data();
        rSAX.push({
            SubDeserializer::onElementBegin,
            SubDeserializer::onText,
            SubDeserializer::onElementEnd,
            pSubInstance
        });

        CIE_BEGIN_EXCEPTION_TRACING
        SubDeserializer::onElementBegin(pSubInstance, rSAX, name, attributes);
        CIE_END_EXCEPTION_TRACING
    } else {
        CIE_THROW(
            Exception,
            "Expecting a \"data\" graphml element on an \"edge\", but got \"" << name << "\"."
        )
    }
}


template <class T>
requires std::is_same_v<typename T::ID,EdgeID>
void GraphML::Deserializer<T>::onText([[maybe_unused]] void*,
                                      [[maybe_unused]] Ref<GraphML::SAXHandler>,
                                      [[maybe_unused]] std::string_view)
{
    CIE_THROW(
        Exception,
        "Found unexpected text data while parsing a node in GraphML."
    )
}


template <class T>
requires std::is_same_v<typename T::ID,EdgeID>
void GraphML::Deserializer<T>::onElementEnd([[maybe_unused]] void* pInstance,
                                            [[maybe_unused]] Ref<GraphML::SAXHandler> rSAX,
                                            [[maybe_unused]] std::string_view name) noexcept
{}


inline void GraphML::Deserializer<void>::onElementBegin([[maybe_unused]] void* pInstance,
                                                        [[maybe_unused]] Ref<GraphML::SAXHandler> rSAX,
                                                        [[maybe_unused]] std::string_view name,
                                                        [[maybe_unused]] std::span<GraphML::AttributePair> attributes) noexcept
{
}


inline void GraphML::Deserializer<void>::onText([[maybe_unused]] void*,
                                                [[maybe_unused]] Ref<GraphML::SAXHandler>,
                                                [[maybe_unused]] std::string_view)
{
    CIE_THROW(
        Exception,
        "Found unexpected text data while parsing void in GraphML."
    )
}


inline void GraphML::Deserializer<void>::onElementEnd([[maybe_unused]] void* pInstance,
                                                      [[maybe_unused]] Ref<GraphML::SAXHandler> rSAX,
                                                      [[maybe_unused]] std::string_view name) noexcept
{
}


inline void GraphML::Deserializer<std::string>::onElementBegin([[maybe_unused]] void* pInstance,
                                                               [[maybe_unused]] Ref<GraphML::SAXHandler> rSAX,
                                                               [[maybe_unused]] std::string_view name,
                                                               [[maybe_unused]] std::span<GraphML::AttributePair> attributes) noexcept
{
}


inline void GraphML::Deserializer<std::string>::onText([[maybe_unused]] void* pInstance,
                                                       [[maybe_unused]] Ref<GraphML::SAXHandler>,
                                                       [[maybe_unused]] std::string_view data)
{
    *static_cast<std::string*>(pInstance) = data;
}


inline void GraphML::Deserializer<std::string>::onElementEnd([[maybe_unused]] void* pInstance,
                                                             [[maybe_unused]] Ref<GraphML::SAXHandler> rSAX,
                                                             [[maybe_unused]] std::string_view name) noexcept
{
}


} // namespace cie::fem::io


#endif
