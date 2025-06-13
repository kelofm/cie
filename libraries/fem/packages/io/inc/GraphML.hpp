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



namespace cie::concepts {
template <class TDeserializer>
concept GraphMLDeserializer
= requires (void* pVoid,
            std::basic_string_view<unsigned char> view,
            std::span<std::pair<std::basic_string_view<unsigned char>,std::basic_string_view<unsigned char>>> views) {
    {TDeserializer::onElementBegin(pVoid, pVoid, view, views)}  -> std::same_as<void>;
    {TDeserializer::onElementEnd(pVoid, pVoid, view)}           -> std::same_as<void>;
    {TDeserializer::onText(pVoid, pVoid, view)}                 -> std::same_as<void>;
};
} // namespace cie::concepts



namespace cie::fem::io {


/// @ingroup fem
class GraphML
{
public:
    using XMLStringView = std::basic_string_view<unsigned char>;

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


    struct DeserializerBase
    {
    public:
        DeserializerBase(Ref<SAXHandler> rSAX);

        virtual ~DeserializerBase() = default;

        Ref<SAXHandler> getSAX() noexcept;

    private:
        Ptr<SAXHandler> _pSAX;
    };


    template <class T>
    struct Deserializer : public GraphML::DeserializerBase {};


    class Input
    {}; // class Input


    class SAXHandler final
    {
    public:
        typedef void (*OnElementBeginCallback)(void*,
                                               void*,
                                               XMLStringView,
                                               std::span<std::pair<XMLStringView,XMLStringView>>);

        typedef void (*OnTextCallback)(void*, void*, XMLStringView);

        typedef void (*OnElementEndCallback)(void*, void*, XMLStringView);

        using State = std::tuple<
            OnElementBeginCallback,
            OnTextCallback,
            OnElementEndCallback,
            void*,
            void*
        >;

        ~SAXHandler();

        void push(State state);

        State pop();

    private:
        SAXHandler(std::istream& rStream);

        SAXHandler(SAXHandler&&) = delete;

        SAXHandler(const SAXHandler&) = delete;

        friend class GraphML::Input;

        struct Impl;
        std::unique_ptr<Impl> _pImpl;
    }; // class SAXHandler


    template <class T>
    struct Serializer
    {
        void header(Ref<XMLElement>)
        requires std::is_same_v<T,std::monostate> {};

        void operator()(Ref<XMLElement>, Ref<const T>)
        requires std::is_same_v<T,std::monostate> {};
    }; // class Serializer


    class Output
    {
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


} // namespace cie::fem::io

#include "packages/io/impl/GraphML_impl.hpp"

#endif
