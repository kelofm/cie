#ifndef CIE_FEM_GRAPHML_HPP
#define CIE_FEM_GRAPHML_HPP

// --- FEM Includes ---
#include "packages/graph/inc/Graph.hpp" // Graph

// --- STL Includes ---
#include <iosfwd> // std::istream
#include <memory> // std::unique_ptr
#include <variant> // std::monostate
#include <filesystem> // std::filesystem::path
#include <span> // std::span


namespace cie::fem::io {


/// @ingroup fem
class GraphML
{
public:
    using XMLStringView = std::span<const unsigned char>;

    using AttributePair = std::pair<std::string_view,std::string_view>;

    class XMLElement final
    {
    public:
        XMLElement(XMLElement&&) noexcept;

        ~XMLElement();

        void addAttribute(std::string_view key,
                          std::string_view value);

        void setValue(std::string_view value);

        XMLElement addChild(std::string_view name);

    private:
        friend class GraphML;

        XMLElement();

        XMLElement(void* pWrapped);

        struct Impl;
        std::unique_ptr<Impl> _pImpl;
    }; // class XMLElement


    class SAXHandler;


    template <class T>
    struct Deserializer {};


    template <class T>
    class DeserializerBase
    {
    public:
        [[nodiscard]] static Ptr<Deserializer<T>> make(typename VoidSafe<T,std::nullptr_t>::Ref rInstance,
                                                       Ref<SAXHandler> rSAX,
                                                       std::string_view elementName);

    protected:
        DeserializerBase(typename VoidSafe<T,std::nullptr_t>::Ref rInstance,
                         Ref<SAXHandler> rSAX) noexcept;

        typename VoidSafe<T,std::nullptr_t>::Ref instance() noexcept;

        Ref<SAXHandler> sax() noexcept;

        template <class TDerived>
        static void release(Ptr<TDerived> pThis,
                            std::string_view elementName) noexcept;

    private:
        Ptr<T> _pInstance;

        Ptr<SAXHandler> _pSAX;
    }; // class DeserializerBase

    class Input
    {
    public:
        Input();

        Input(Input&&) noexcept;

        Input(Ref<std::istream> rStream);

        ~Input();

        template <class TVertexData, class TEdgeData, class TGraphData>
        void operator()(Ref<Graph<TVertexData,TEdgeData,TGraphData>> rGraph);

    private:
        Ref<std::istream> stream() noexcept;

        struct Impl;
        std::unique_ptr<Impl> _pImpl;
    }; // class Input


    class SAXHandler final
    {
    public:
        typedef void (*OnElementBeginCallback)(void*,
                                               std::string_view,
                                               std::span<AttributePair>);

        typedef void (*OnTextCallback)(void*, std::string_view);

        typedef void (*OnElementEndCallback)(void*, std::string_view);

        using State = std::tuple<
            void*,
            OnElementBeginCallback,
            OnTextCallback,
            OnElementEndCallback
        >;

        ~SAXHandler();

        void push(State state);

    private:
        SAXHandler(std::istream& rStream);

        SAXHandler(SAXHandler&&) = delete;

        SAXHandler(const SAXHandler&) = delete;

        void parse(std::size_t bufferSize = 0x800000);

        friend class GraphML::Input;

        struct Impl;
        std::unique_ptr<Impl> _pImpl;
    }; // class SAXHandler


    template <class T>
    struct Serializer {
        void header(Ref<XMLElement>)
        requires std::is_same_v<T,std::monostate> {};

        void operator()(Ref<XMLElement>, Ref<const T>)
        requires std::is_same_v<T,std::monostate> {};
    }; // class Serializer


    class Output {
    public:
        Output();

        Output(Ref<const std::filesystem::path> rOutputPath);

        ~Output();

        template <class TVertexData, class TEdgeData, class TGraphData>
        void operator()(Ref<const Graph<TVertexData,TEdgeData,TGraphData>> rGraph);

    private:
        XMLElement root();

        void write();

        struct Impl;
        std::unique_ptr<Impl> _pImpl;
    }; // class Output
}; // class GraphML


template <class TVertexData, class TEdgeData, class TGraphData>
struct GraphML::Serializer<Graph<TVertexData,TEdgeData,TGraphData>>
{
    void header(Ref<GraphML::XMLElement> rElement) const;

    void operator()(Ref<GraphML::XMLElement> rElement,
                    Ref<const Graph<TVertexData,TEdgeData,TGraphData>> rInstance) const;
}; // struct Serializer<Graph>


template <class TVertexData, class TEdgeData, class TGraphData>
struct GraphML::Deserializer<Graph<TVertexData,TEdgeData,TGraphData>>
    : public GraphML::DeserializerBase<Graph<TVertexData,TEdgeData,TGraphData>>
{
    using GraphML::DeserializerBase<Graph<TVertexData,TEdgeData,TGraphData>>::DeserializerBase;

    static void onElementBegin(void* pInstance,
                               std::string_view name,
                               std::span<GraphML::AttributePair> attributes);

    static void onText(void* pInstance,
                       std::string_view data);

    static void onElementEnd(void* pInstance,
                             std::string_view name);
}; // struct Deserializer<Graph>


template <class T>
requires std::is_same_v<typename T::ID,VertexID>
struct GraphML::Deserializer<T>
    : public GraphML::DeserializerBase<T>
{
    using GraphML::DeserializerBase<T>::DeserializerBase;

    static void onElementBegin(Ptr<void> pThis,
                               std::string_view name,
                               std::span<GraphML::AttributePair> attributes);

    static void onText(Ptr<void> pThis,
                       std::string_view data);

    static void onElementEnd(Ptr<void> pThis,
                             std::string_view name) noexcept;
}; // struct Deserializer<Graph::Vertex>


template <class T>
requires std::is_same_v<typename T::ID,EdgeID>
struct GraphML::Deserializer<T>
    : public GraphML::DeserializerBase<T>
{
    using GraphML::DeserializerBase<T>::DeserializerBase;

    static void onElementBegin(Ptr<void> pThis,
                               std::string_view name,
                               std::span<GraphML::AttributePair> attributes);

    static void onText(Ptr<void> pThis,
                       std::string_view data);

    static void onElementEnd(Ptr<void> pThis,
                             std::string_view name) noexcept;
}; // struct Deserializer<Graph::Edge>


} // namespace cie::fem::io

#include "packages/io/impl/GraphML_impl.hpp"

#endif
