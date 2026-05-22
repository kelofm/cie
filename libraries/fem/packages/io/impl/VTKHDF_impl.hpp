#pragma once

// --- FEM Includes ---
#include "packages/io/inc/VTKHDF.hpp"
#include "packages/utilities/inc/kernel.hpp"

// --- Utility Includes ---
#include "packages/macros/inc/exceptions.hpp"
#include "packages/maths/inc/power.hpp"
#include "packages/maths/inc/OuterProduct.hpp"
#include "packages/concurrency/inc/ParallelFor.hpp"


namespace cie::io {


template <class TValue>
void VTKHDF::Output::writeFieldVariables(
    std::string_view groupName,
    Ref<const std::vector<std::tuple<
        std::string,
        std::vector<std::size_t>,
        std::span<const TValue>
    >>> fieldVariables) {
        CIE_BEGIN_EXCEPTION_TRACING
            const Prefix groupPrefix = Prefix("/VTKHDF") / groupName / "PointData";
            this->makeGroup(groupPrefix);
            for (const auto& [rName, rShape, rValues] : fieldVariables) {
                switch (rShape.size()) {
                    case 0:
                        this->writeDataset<TValue,1>(
                            groupPrefix / rName,
                            {rValues.size()},
                            rValues);
                        break;
                    case 1:
                        this->writeDataset<TValue,2>(
                            groupPrefix / rName,
                            {rValues.size() / rShape.front(), rShape.front()},
                            rValues);
                        break;
                    default: CIE_THROW(Exception, "unhandled field variable array rank " << rShape.size())
                }
            } // for rName, rShape, rValues in fieldVariables
        CIE_END_EXCEPTION_TRACING
}


template <class TValue>
void VTKHDF::Output::writeCellVariables(
    std::string_view groupName,
    std::vector<std::tuple<
        std::string,
        std::vector<std::size_t>,
        std::span<const TValue>
    >> cellVariables) {
        CIE_BEGIN_EXCEPTION_TRACING
            const Prefix groupPrefix = Prefix("/VTKHDF") / groupName / "CellData";
            CIE_BEGIN_EXCEPTION_TRACING
                this->makeGroup(groupPrefix);
            CIE_END_EXCEPTION_TRACING

            for (const auto& [rName, rShape, rValues] : cellVariables) {
                const std::size_t componentCount = std::accumulate(
                    rShape.begin(),
                    rShape.end(),
                    1ul,
                    std::multiplies<std::size_t>());
                const std::size_t cellCount = rValues.size() / componentCount;
                CIE_CHECK(
                    cellCount * componentCount == rValues.size(),
                    std::format(
                        "expecting cell array \"{}\" to be of size {}, but has {} components",
                        rName, cellCount * componentCount, rValues.size()))

                switch (rShape.size()) {
                    case 0:
                        this->writeDataset<TValue,1>(
                            groupPrefix / rName,
                            {cellCount},
                            rValues);
                        break;
                    case 1:
                        this->writeDataset<TValue,2>(
                            groupPrefix / rName,
                            {cellCount, componentCount},
                            rValues);
                        break;
                    default: CIE_THROW(Exception, "unexpected cell variable array rank " << rShape.size())
                }; // switch rShape.size()
            } // for rName, rShape, rValues in cellVariables
        CIE_END_EXCEPTION_TRACING
}


template <fem::CellLike TCell>
void VTKHDF::Output::writeCells(
    std::string_view groupName,
    std::span<const TCell> cells) {
        using Scalar = typename TCell::Value;
        mp::ThreadPoolBase threads;
        const Prefix meshPrefix = Prefix("/VTKHDF") / groupName;

        CIE_BEGIN_EXCEPTION_TRACING
            this->makeGroup(meshPrefix);
            this->writeAttribute(
                meshPrefix,
                "Type",
                "UnstructuredGrid");
            const std::array<int,2> version {2, 4};
            this->writeAttribute<int>(
                meshPrefix,
                "Version",
                version);
        CIE_END_EXCEPTION_TRACING

        CIE_BEGIN_EXCEPTION_TRACING
            // Dataset sizes.
            std::size_t cellCount = cells.size();
            this->writeDataset<std::size_t,1>(
                meshPrefix / "NumberOfCells",
                {1},
                {&cellCount, 1});

            constexpr std::size_t pointsPerCell = intPow<std::size_t>(2, TCell::ParametricDimension);
            const std::size_t topologySize = cellCount * pointsPerCell;
            this->writeDataset<std::size_t,1>(
                meshPrefix / "NumberOfConnectivityIds",
                {1},
                {&topologySize, 1});

            const std::size_t pointCount = pointsPerCell * cellCount;
            this->writeDataset<std::size_t,1>(
                meshPrefix / "NumberOfPoints",
                {1},
                {&pointCount, 1});

            // Cell-wise datasets.
            CIE_BEGIN_EXCEPTION_TRACING
                std::vector<int> topology(topologySize);
                std::vector<std::uint8_t> topologyTypes(cellCount);
                std::vector<int> offsets(cellCount + 1);
                std::vector<Scalar> points(pointCount * 3);
                offsets.front() = 0ul;

                const auto kernel = [&] (
                    std::size_t iCell,
                    Ref<std::vector<typename TCell::Value>> rBuffer) {
                        Ref<const TCell> rCell = cells[iCell];
                        const std::size_t iCellBegin = pointsPerCell * iCell;

                        // Get topology.
                        if constexpr (TCell::ParametricDimension == 1) {
                            topology[iCellBegin    ] = iCellBegin    ;
                            topology[iCellBegin + 1] = iCellBegin + 1;
                        } else if constexpr (TCell::ParametricDimension == 2) {
                            topology[iCellBegin    ] = iCellBegin    ;
                            topology[iCellBegin + 1] = iCellBegin + 1;
                            topology[iCellBegin + 2] = iCellBegin + 3;
                            topology[iCellBegin + 3] = iCellBegin + 2;
                        } else if constexpr (TCell::ParametricDimension == 3) {
                            topology[iCellBegin    ] = iCellBegin    ;
                            topology[iCellBegin + 1] = iCellBegin + 1;
                            topology[iCellBegin + 2] = iCellBegin + 3;
                            topology[iCellBegin + 3] = iCellBegin + 2;
                            topology[iCellBegin + 4] = iCellBegin + 4;
                            topology[iCellBegin + 5] = iCellBegin + 5;
                            topology[iCellBegin + 6] = iCellBegin + 7;
                            topology[iCellBegin + 7] = iCellBegin + 6;
                        } else static_assert(TCell::ParametricDimension == 1, "unsupported dimension");

                        // Get topology type.
                        if constexpr (TCell::ParametricDimension == 1) {
                            topologyTypes[iCell] = /*VTK_LINE*/ 3;
                        } else if constexpr (TCell::ParametricDimension == 2) {
                            topologyTypes[iCell] = /*VTK_QUAD*/ 9;
                        } else if constexpr (TCell::ParametricDimension == 3) {
                            topologyTypes[iCell] = /*VTK_HEXAHEDRON*/ 12;
                        } else static_assert(TCell::ParametricDimension == 1, "unsupported dimension");

                        // Get the index of the cell's topology begin.
                        offsets[iCell + 1] = (iCell + 1) * pointsPerCell;

                        // Get node coordinates.
                        std::array<fem::ParametricCoordinate<Scalar>,TCell::ParametricDimension> parametricCoordinates;
                        std::array<fem::PhysicalCoordinate<Scalar>,TCell::PhysicalDimension> physicalCoordinates;
                        std::array<std::uint8_t,TCell::ParametricDimension> state;
                        std::fill_n(
                            state.data(),
                            state.size(),
                            static_cast<std::uint8_t>(0));

                        std::size_t iComponentBegin = iCell * pointsPerCell * 3;
                        rBuffer.resize(rCell.makeSpatialTransform().bufferSize());
                        do {
                            std::transform(
                                state.begin(),
                                state.end(),
                                parametricCoordinates.begin(),
                                [] (std::uint8_t state) -> fem::ParametricCoordinate<Scalar> {
                                    return state ? static_cast<Scalar>(1) : static_cast<Scalar>(-1);
                                });

                            rCell.transform(
                                fem::Kernel<TCell::ParametricDimension,Scalar>::view(parametricCoordinates),
                                fem::Kernel<TCell::PhysicalDimension,Scalar>::view(physicalCoordinates),
                                rBuffer);

                            std::copy_n(
                                physicalCoordinates.data(),
                                std::min<std::size_t>(3, TCell::PhysicalDimension),
                                points.data() + iComponentBegin);
                            std::fill_n(
                                points.data() + iComponentBegin + std::min<std::size_t>(3, TCell::PhysicalDimension),
                                3 - std::min<std::size_t>(3, TCell::PhysicalDimension),
                                static_cast<Scalar>(0));

                            iComponentBegin += 3;
                        } while (cie::maths::OuterProduct<TCell::ParametricDimension>::next(2u, state.data()));
                }; // kernel

                mp::ParallelFor<>(threads).firstPrivate(std::vector<typename TCell::Value>()).execute(
                    cells.size(),
                    [&] (std::size_t iCell, Ref<std::vector<typename TCell::Value>> rBuffer) {
                        kernel(iCell, rBuffer);
                    }); // for rCell in rCells


                CIE_BEGIN_EXCEPTION_TRACING
                    this->writeDataset<int,1>(
                        meshPrefix / "Connectivity",
                        {topology.size()},
                        topology);
                CIE_END_EXCEPTION_TRACING

                CIE_BEGIN_EXCEPTION_TRACING
                    this->writeDataset<std::uint8_t,1>(
                        meshPrefix / "Types",
                        {topologyTypes.size()},
                        topologyTypes);
                CIE_END_EXCEPTION_TRACING

                CIE_BEGIN_EXCEPTION_TRACING
                    this->writeDataset<int,1>(
                        meshPrefix / "Offsets",
                        {offsets.size()},
                        offsets);
                CIE_END_EXCEPTION_TRACING

                CIE_BEGIN_EXCEPTION_TRACING
                    this->writeDataset<Scalar,2>(
                        meshPrefix / "Points",
                        {pointCount, 3},
                        points);
                CIE_END_EXCEPTION_TRACING
        CIE_END_EXCEPTION_TRACING

        CIE_END_EXCEPTION_TRACING

        CIE_BEGIN_EXCEPTION_TRACING
            this->makeGroup(Prefix("/VTKHDF") / "Assembly" / groupName);
            this->link(
                Prefix("/VTKHDF") / "Assembly" / groupName / groupName,
                meshPrefix);
        CIE_END_EXCEPTION_TRACING
}


template <class TValue, fem::DiscretizationLike TMesh, fem::CellLike TCell>
requires std::is_same_v<typename TCell::Value,TValue>
void VTKHDF::Output::writeFieldVariables(
    std::string_view groupName,
    Ref<const TMesh> rMesh,
    std::span<const TCell> cells,
    Ref<const fem::Assembler> rAssembler,
    std::vector<std::pair<
        std::string,
        std::span<const TValue>>
    > fieldVariables) {
        if (fieldVariables.empty()) return;

        const Prefix groupPrefix = Prefix("/VTKHDF") / groupName / "PointData";
        const std::size_t dofCount = rAssembler.dofCount();

        // Sanity checks.
        CIE_CHECK(
            cells.size() == rMesh.vertices().size(),
            std::format(
                "number of cells in the mesh({}) does not match the number of provided cells ({})",
                rMesh.vertices().size(),
                cells.size()))
        const std::size_t cellCount = cells.size();
        constexpr std::size_t pointsPerCell = intPow<std::size_t>(2, TCell::ParametricDimension);
        const std::size_t pointCount = cellCount * pointsPerCell;
        std::vector<std::size_t> componentCounts(fieldVariables.size());
        std::vector<std::vector<TValue>> fieldVariableOutput(fieldVariables.size());

        for (std::size_t iField=0ul; iField<fieldVariables.size(); ++iField) {
            const auto& [rVariableName, rValues] = fieldVariables[iField];
            CIE_CHECK(
                !(rValues.size() % dofCount),
                std::format(
                    "expecting field variable array \"{}\" size to be a multiple of {}, but has {} entries",
                    rVariableName, dofCount, rValues.size()))

            const std::size_t componentCount = rValues.size() / dofCount;
            componentCounts[iField] = componentCount;
            fieldVariableOutput[iField].resize(pointCount * componentCount);
        } // for iField in fieldVariables.size()

        CIE_BEGIN_EXCEPTION_TRACING
            this->makeGroup(groupPrefix);
        CIE_END_EXCEPTION_TRACING

        CIE_BEGIN_EXCEPTION_TRACING
            const auto kernel = [&] (
                std::size_t iCell,
                Ref<std::vector<TValue>> rAnsatzBuffer,
                Ref<std::vector<TValue>> rExpressionBuffer) -> void {
                    Ref<const TCell> rCell = cells[iCell];
                    const auto& rGlobalDofIndices = rAssembler[rCell.id()];
                    const auto& rAnsatz = rMesh.data().ansatz(rCell.ansatzID());
                    rAnsatzBuffer.resize(rAnsatz.size());
                    rExpressionBuffer.resize(rAnsatz.bufferSize());

                    std::array<TValue,TCell::ParametricDimension> parametricCoordinates;
                    std::array<std::uint8_t,TCell::ParametricDimension> state;
                    std::fill_n(
                        state.data(),
                        state.size(),
                        static_cast<std::uint8_t>(0));

                    std::size_t iLocalPoint = 0ul;
                    do {
                        // Construct the node's parametric coordinates.
                        std::transform(
                            state.begin(),
                            state.end(),
                            parametricCoordinates.begin(),
                            [] (std::uint8_t state) -> TValue {
                                return state ? static_cast<TValue>(1) : static_cast<TValue>(-1);
                            });

                        // Evaluate the ansatz space.
                        rAnsatz.evaluate(
                            parametricCoordinates,
                            rAnsatzBuffer,
                            rExpressionBuffer);

                        // Compute field values.
                        const std::size_t iPoint = iCell * pointsPerCell + iLocalPoint;
                        for (std::size_t iField=0ul; iField<fieldVariables.size(); ++iField) {
                            const std::size_t iOutputComponentBegin = iPoint * componentCounts[iField];
                            for (std::size_t iComponent=0ul; iComponent<componentCounts[iField]; ++iComponent) {
                                Ref<TValue> rOutputComponent = fieldVariableOutput[iField][iOutputComponentBegin + iComponent];
                                rOutputComponent = static_cast<TValue>(0);
                                for (std::size_t iAnsatz=0ul; iAnsatz<rAnsatzBuffer.size(); ++iAnsatz) {
                                    const std::size_t iGlobalDof = componentCounts[iField] * rGlobalDofIndices[iAnsatz] + iComponent;
                                    rOutputComponent += fieldVariables[iField].second[iGlobalDof] * rAnsatzBuffer[iAnsatz];
                                } // for iAnsatz in range(rAnsatzBuffer.size())
                            } // for iComponent in range(componentCounts.size())
                        } // for iField in range(fieldVariables.size())

                        ++iLocalPoint;
                    } while (cie::maths::OuterProduct<TCell::ParametricDimension>::next(2u, state.data()));
            }; // kernel

            // Run the kernel on the cells.
            mp::ThreadPoolBase threads;
            mp::ParallelFor<>(threads).firstPrivate(std::vector<TValue>(), std::vector<TValue>()).execute(
                cells.size(),
                [&kernel] (std::size_t iCell, Ref<std::vector<TValue>> rAnsatzBuffer, Ref<std::vector<TValue>> rExpressionBuffer) {
                    kernel(iCell, rAnsatzBuffer, rExpressionBuffer);});

            // Write fields to the HDF5 file as PointData.
            for (std::size_t iField=0ul; iField<fieldVariables.size(); ++iField) {
                const std::size_t componentCount = componentCounts[iField];
                std::string_view fieldName = fieldVariables[iField].first;
                switch (componentCount) {
                    case 1:
                        this->writeDataset<TValue,1>(
                            groupPrefix / fieldName,
                            {pointCount},
                            fieldVariableOutput[iField]);
                        break;
                    case 2:
                        this->writeDataset<TValue,2>(
                            groupPrefix / fieldName,
                            {pointCount, componentCount},
                            fieldVariableOutput[iField]);
                        break;
                    default:
                        CIE_THROW(Exception, "unexpected field array rank " << componentCount)
                }
            } // for iField in range(fieldVariables.size())
        CIE_END_EXCEPTION_TRACING
}


} // namespace cie::io
