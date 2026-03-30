#pragma once

#ifdef CIE_ENABLE_HDF5

// --- FEM Includes ---
#include "packages/numeric/inc/CellBase.hpp"
#include "packages/numeric/inc/MeshBase.hpp"
#include "packages/graph/inc/Assembler.hpp"

// --- STL Includes ---
#include <filesystem> // std::filesystem::path
#include <memory> // std::unique_ptr
#include <string_view> // std::string_view
#include <optional> // std::optional
#include <vector> // std::vector


namespace cie::io {


struct VTKHDF {
    class Output {
    public:
        enum class AttributeType {
            Mesh,
            Cell,
            Vertex
        }; // AttributeType

        Output();

        Output(Ref<const std::filesystem::path> rPath);

        ~Output();

        /// @brief Write a mesh of cells to an unstructured grid.
        template <fem::CellLike TCell>
        void operator()(
            std::string_view groupName,
            std::span<const TCell> cells);

        /// @brief Write a mesh of cells to an unstructured grid.
        template <fem::DiscretizationLike TMesh>
        void operator()(
            std::string_view groupName,
            Ref<const TMesh> rMesh);

        template <class TValue, unsigned Dimension>
        void writePointCloud(
            std::string_view groupName,
            std::span<const TValue> coordinates,
            std::size_t gridSize = 0ul);

        template <class TValue, fem::DiscretizationLike TMesh>
        requires std::is_same_v<typename TMesh::Vertex::Data::Value,TValue>
        void writeFieldVariables(
            std::string_view groupName,
            Ref<const TMesh> rMesh,
            Ref<const fem::Assembler> rAssembler,
            std::vector<std::pair<
                std::string,
                std::span<const TValue>>
            > fieldVariables);

        template <class TValue, fem::DiscretizationLike TMesh, fem::CellLike TCell>
        requires (
            std::is_same_v<typename TMesh::Vertex::Data,TCell>
            && std::is_same_v<typename TCell::Value,TValue>)
        void writeFieldVariables(
            std::string_view groupName,
            Ref<const TMesh> rMesh,
            std::span<const TCell> cells,
            Ref<const fem::Assembler> rAssembler,
            std::vector<std::pair<
                std::string,
                std::span<const TValue>>
            > fieldVariables);

        template <class TValue>
        void writeFieldVariables(
            std::string_view groupName,
            Ref<const std::vector<std::tuple<
                std::string,
                std::vector<std::size_t>,
                std::span<const TValue>
            >>> fieldVariables);

        template <class TValue>
        void writeCellVariables(
            std::string_view groupName,
            std::vector<std::tuple<
                std::string,
                std::vector<std::size_t>,
                std::span<const TValue>
            >> cellVariables);

    protected:
        using Prefix = std::filesystem::path;

        template <class TMesh, class TCell>
        void writeMesh(
            std::optional<std::reference_wrapper<const TMesh>> rMaybeMesh,
            std::optional<std::span<const TCell>> maybeCells,
            std::string_view groupName);

        template <class TMesh, class TCell, class TValue>
        void writeFieldVariablesImpl(
            std::string_view groupName,
            Ref<const TMesh> rMesh,
            std::optional<std::span<const TCell>> maybeCells,
            Ref<const fem::Assembler> rAssembler,
            std::vector<std::pair<
                std::string,
                std::span<const TValue>>
            > fieldVariables);

        void makeGroup(Ref<const Prefix> rPrefix);

        void link(
            Ref<const Prefix> rSource,
            Ref<const Prefix> rTarget);

        void writeAttribute(
            Ref<const Prefix> rPrefix,
            Ref<const std::string> rName,
            std::string_view value);

        template <class T>
        void writeAttribute(
            Ref<const Prefix> rPrefix,
            Ref<const std::string> rName,
            std::span<const T> value);

        template <class TValue, unsigned Dimension>
        void writeDataset(
            Ref<const Prefix> rPrefix,
            std::array<const std::size_t,Dimension> shape,
            std::span<const TValue> data);

    private:
        struct Impl;
        std::unique_ptr<Impl> _pImpl;
    }; // class Output
}; // struct VTKHDF


} // namespace cie::io

#include "packages/io/impl/VTKHDF_impl.hpp"

#endif
