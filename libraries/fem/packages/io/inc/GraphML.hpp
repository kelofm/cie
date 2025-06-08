#ifndef CIE_FEM_GRAPHML_HPP
#define CIE_FEM_GRAPHML_HPP

// --- FEM Includes ---
#include "packages/graph/inc/Graph.hpp" // Graph

// --- STL Includes ---
#include <string> // std::string
#include <memory> // std::unique_ptr
#include <variant> // std::monostate
#include <filesystem> // std::filesystem::path


namespace cie::fem::io {


/// @ingroup fem
class GraphML
{
public:
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

    private:
        struct Impl;
        std::unique_ptr<Impl> _pImpl;
    }; // class XMLElement


    template <class T>
    struct Serializer
    {
        void header(Ref<XMLElement>)
        requires std::is_same_v<T,std::monostate> {};

        void operator()(Ref<XMLElement>, Ref<const T>)
        requires std::is_same_v<T,std::monostate> {};
    }; // class Serializer


    class Input
    {}; // class Input


    class Output
    {
    public:
        Output();

        Output(Ref<const std::filesystem::path> rOutputPath);

        ~Output();

        template <class TVertexData, class TEdgeData>
        void operator()(Ref<const Graph<TVertexData,TEdgeData>> rGraph);

    private:
        XMLElement root();

        void write();

    private:
        struct Impl;
        std::unique_ptr<Impl> _pImpl;
    }; // class Output
}; // class GraphML


} // namespace cie::fem::io

#include "packages/io/impl/GraphML_impl.hpp"

#endif
