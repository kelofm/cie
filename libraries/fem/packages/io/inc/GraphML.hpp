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


    template <class T>
    struct Deserializer {};


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
                                               Ref<SAXHandler>,
                                               XMLStringView,
                                               std::span<std::pair<XMLStringView,XMLStringView>>);

        typedef void (*OnTextCallback)(void*,
                                       Ref<SAXHandler>,
                                       XMLStringView);

        typedef void (*OnElementEndCallback)(void*,
                                             Ref<SAXHandler>,
                                             XMLStringView);

        using State = std::tuple<
            OnElementBeginCallback,
            OnTextCallback,
            OnElementEndCallback,
            void*
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
