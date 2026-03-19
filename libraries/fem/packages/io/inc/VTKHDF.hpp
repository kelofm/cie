#pragma once

#ifdef CIE_ENABLE_HDF5

// --- FEM Includes ---
#include "packages/numeric/inc/CellBase.hpp"
#include "packages/numeric/inc/MeshBase.hpp"

// --- STL Includes ---
#include <filesystem> // std::filesystem::path
#include <memory> // std::unique_ptr


namespace cie::io {


struct VTKHDF {
    class Output {
    public:
        Output();

        Output(Ref<const std::filesystem::path> rPath);

        ~Output();

        template <fem::DiscretizationLike TMesh>
        void operator()(Ref<const TMesh>);

    private:
        using Prefix = std::filesystem::path;

        void writeAttribute(
            Ref<const Prefix> rPrefix,
            Ref<const std::string> rName,
            std::string_view value);

        template <class T>
        void writeAttribute(
            Ref<const Prefix> rPrefix,
            Ref<const std::string> rName,
            std::span<const T> value);

        struct Impl;
        std::unique_ptr<Impl> _pImpl;
    }; // class Output
}; // struct VTKHDF


} // namespace cie::io

#include "packages/io/impl/VTKHDF_impl.hpp"

#endif
