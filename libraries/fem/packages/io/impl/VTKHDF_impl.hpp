#pragma once

// --- FEM Includes ---
#include "packages/io/inc/VTKHDF.hpp"
#include "packages/utilities/inc/kernel.hpp"

// --- Utility Includes ---
#include "packages/macros/inc/exceptions.hpp"
#include "packages/maths/inc/power.hpp"
#include "packages/maths/inc/OuterProduct.hpp"


namespace cie::io {


template <fem::DiscretizationLike TMesh>
void VTKHDF::Output::operator()(
    Ref<const TMesh> rMesh,
    std::string_view name) {
        using Scalar = typename TMesh::Vertex::Data::Value;
        constexpr unsigned Dimension = TMesh::Data::Ansatz::Dimension;

        const Prefix meshPrefix = Prefix("/VTKHDF") / name;
        const auto& rCells = rMesh.vertices();

        CIE_BEGIN_EXCEPTION_TRACING
            this->makeGroup(meshPrefix);
            this->makeGroup(Prefix("/VTKHDF/Assembly") / name);
            this->link(
                Prefix("/VTKHDF/Assembly") / name / name,
                meshPrefix);
            this->writeAttribute(
                meshPrefix,
                "Type",
                "UnstructuredGrid");
        CIE_END_EXCEPTION_TRACING

        CIE_BEGIN_EXCEPTION_TRACING
            // Dataset sizes.
            const std::size_t cellCount = rCells.size();
            this->writeDataset<std::size_t,1>(
                meshPrefix / "NumberOfCells",
                {1},
                {&cellCount, 1});

            constexpr std::size_t pointsPerCell = intPow<std::size_t>(2, Dimension);
            this->writeDataset<std::size_t,1>(
                meshPrefix / "NumberOfConnectivityIds",
                {1},
                {&cellCount, 1});

            const std::size_t pointCount = pointsPerCell * cellCount;
            this->writeDataset<std::size_t,1>(
                meshPrefix / "NumberOfPoints",
                {1},
                {&pointCount, 1});

            // Cell-wise datasets.
            CIE_BEGIN_EXCEPTION_TRACING

                std::vector<std::size_t> topology(cellCount * pointsPerCell);
                std::vector<Scalar> points(pointCount * 3);
                std::size_t iCell = 0ul;

                for (Ref<const typename TMesh::Vertex> rCell : rCells) {
                        const std::size_t iCellBegin = pointsPerCell * iCell++;

                        if constexpr (Dimension == 1) {
                            topology[iCellBegin] = iCellBegin;
                            topology[iCellBegin + 1] = iCellBegin + 1;
                        } else if constexpr (Dimension == 2) {
                            topology[iCellBegin] = iCellBegin;
                            topology[iCellBegin + 1] = iCellBegin + 1;
                            topology[iCellBegin + 2] = iCellBegin + 3;
                            topology[iCellBegin + 3] = iCellBegin + 2;
                        } else if constexpr (Dimension == 3) {
                            topology[iCellBegin] = iCellBegin;
                            topology[iCellBegin + 1] = iCellBegin + 1;
                            topology[iCellBegin + 2] = iCellBegin + 3;
                            topology[iCellBegin + 3] = iCellBegin + 2;
                            topology[iCellBegin + 4] = iCellBegin + 4;
                            topology[iCellBegin + 5] = iCellBegin + 5;
                            topology[iCellBegin + 6] = iCellBegin + 7;
                            topology[iCellBegin + 7] = iCellBegin + 6;
                        } else static_assert(Dimension == 1, "unsupported dimension");

                        std::array<fem::ParametricCoordinate<Scalar>,TMesh::Vertex::Data::ParametricDimension> parametricCoordinates;
                        std::array<fem::PhysicalCoordinate<Scalar>,TMesh::Vertex::Data::PhysicalDimension> physicalCoordinates;
                        std::array<std::uint8_t,intPow(2,TMesh::Vertex::Data::ParametricDimension)> state;
                        std::fill_n(
                            state.data(),
                            state.size(),
                            static_cast<std::uint8_t>(0));

                        std::size_t iComponentBegin = iCell * pointsPerCell * 3;
                        do {
                            std::transform(
                                state.begin(),
                                state.end(),
                                parametricCoordinates.begin(),
                                [] (std::uint8_t state) -> fem::ParametricCoordinate<Scalar> {
                                    return state ? static_cast<Scalar>(1) : static_cast<Scalar>(-1);
                                });

                            rCell.data().transform(
                                fem::Kernel<TMesh::Vertex::Data::ParametricDimension,Scalar>::view(parametricCoordinates),
                                fem::Kernel<TMesh::Vertex::Data::PhysicalDimension,Scalar>::view(physicalCoordinates));

                            std::copy_n(
                                physicalCoordinates.data(),
                                3,
                                points.data() + iComponentBegin);
                            iComponentBegin += 3;
                        } while (cie::maths::OuterProduct<TMesh::Vertex::Data::ParametricDimension>::next(2u, state.data()));
                    } // for rCell in rCells

                    this->writeDataset<std::size_t,1>(
                        meshPrefix / "Connectivity",
                        {topology.size()},
                        topology);

                    this->writeDataset<Scalar,2>(
                        meshPrefix / "Points",
                        {pointCount, 3},
                        points);
            CIE_END_EXCEPTION_TRACING

        CIE_END_EXCEPTION_TRACING
}


} // namespace cie::io
