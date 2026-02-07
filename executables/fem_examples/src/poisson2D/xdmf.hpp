#pragma once

// --- External Includes ---
#ifdef CIE_ENABLE_HDF5
    #include "H5Cpp.h"
#endif

// --- Internal Includes ---
#include "poisson2D/CellData.hpp"
#include "poisson2D/BoundaryData.hpp"
#include "poisson2D/MeshData.hpp"
#include "poisson2D/mesh.hpp"
#include "poisson2D/integration.hpp"
#include "poisson2D/constraints.hpp"

// --- Utility Includes ---
#include "packages/commandline/inc/ArgParse.hpp"

// --- STL Includes ---
#include <filesystem>
#include <ranges> // ranges::iota


namespace cie::fem {


struct Sample {
    StaticArray<CellData::PhysicalCoordinate,Dimension> position;
    Scalar                                              state;
    Scalar                                              load;
    Scalar                                              residual;
    unsigned                                            cellID;
}; // struct Sample


class Output {
public:
    Output(Ref<const std::filesystem::path> rPath = std::filesystem::path())
        : _xdmfFile(rPath / "poisson2D.xdmf") {
            #ifdef CIE_ENABLE_HDF5
                const auto path = (rPath / "poisson2D.h5").string();
                H5std_string h5Path;
                h5Path.resize(path.size());
                std::copy(
                    path.begin(),
                    path.end(),
                    h5Path.begin());
                _h5File = H5::H5File(
                    h5Path,
                    H5F_ACC_TRUNC);
            #endif
        }

    template <class T>
    void write(
        std::span<const T> data,
        [[maybe_unused]] Ref<const std::string> rDatasetName,
        [[maybe_unused]] Ref<const std::vector<int>> rShape) {
        #ifdef CIE_ENABLE_HDF5
            //for (const auto& c : data) _xdmfFile << c << " ";

            H5::PredType valueType = H5::PredType::IEEE_F32LE;
            if constexpr (std::is_same_v<T,int>) {
                if constexpr (sizeof(int) == 4) valueType = H5::PredType::NATIVE_INT32;
                else if constexpr (sizeof(int) == 8) valueType = H5::PredType::NATIVE_INT64;
                else static_assert(std::is_same_v<T,void>, "unsupported type");
            } else if constexpr (std::is_same_v<T,unsigned>) {
                if constexpr (sizeof(unsigned) == 4) valueType = H5::PredType::NATIVE_UINT32;
                else if constexpr (sizeof(unsigned) == 8) valueType = H5::PredType::NATIVE_UINT64;
                else static_assert(std::is_same_v<T,void>, "unsupported type");
            } else if constexpr (std::is_same_v<T,std::size_t>) {
                if constexpr (sizeof(std::size_t) == 4) valueType = H5::PredType::NATIVE_UINT32;
                else if constexpr (sizeof(std::size_t) == 8) valueType = H5::PredType::NATIVE_UINT64;
                else static_assert(std::is_same_v<T,void>, "unsupported type");
            } else if constexpr (std::is_same_v<T,std::int8_t>) valueType = H5::PredType::NATIVE_INT8;
            else if constexpr (std::is_same_v<T,float>) valueType = H5::PredType::IEEE_F32LE;
            else if constexpr (std::is_same_v<T,double>) valueType = H5::PredType::IEEE_F64LE;
            else static_assert(std::is_same_v<T,void>, "unsupported_type");

            std::vector<hsize_t> shape(rShape.begin(), rShape.end());
            H5::DataSpace dataSpace(
                shape.size(),
                shape.data());
            H5::DataSet dataset = _h5File.createDataSet(
                rDatasetName.c_str(),
                valueType,
                dataSpace);
            dataset.write(data.data(), valueType);

            _xdmfFile << "poisson2D.h5:/" << rDatasetName;
        #else
            for (const auto& c : data) rOut._xdmfFile << c << " ";
        #endif
    }

    template <class T>
    friend Ref<Output> operator<<(Ref<Output> rOut, Ref<const T> rMessage) {
        rOut._xdmfFile << rMessage;
        return rOut;
    }

private:
    std::ofstream _xdmfFile;

    #ifdef CIE_ENABLE_HDF5
        H5::H5File _h5File;
    #endif
}; // class Output


void postprocess(
    CSRWrapper lhs,
    std::span<const Scalar> solution,
    std::span<const Scalar> rhs,
    Ref<const Mesh> rMesh,
    std::span<const CellData> contiguousCellData,
    std::span<const BoundarySegment> boundarySegments,
    Ref<const BVH> rBVH,
    Ref<const Assembler> rAssembler,
    Ref<const utils::ArgParse::Results> rArguments,
    Ref<mp::ThreadPoolBase> rThreads)
{
    CIE_BEGIN_EXCEPTION_TRACING
    const unsigned postprocessResolution = rArguments.get<std::size_t>("scatter-resolution");
    Output output;

    // Collect output for the bounding volume hierarchy.
    DynamicArray<StaticArray<Scalar,Dimension*intPow(2,Dimension)>> boundingVolumes; // {p0x, p0y, p1x, p1y, p2x, p2y, p3x, p3y}
    rBVH.visit([&boundingVolumes](const auto& rBox) -> bool {
        constexpr unsigned cornerCount = intPow(2,Dimension);
        boundingVolumes.emplace_back();
        geo::Box<Dimension,Scalar>::makeCorners(
            rBox.base(),
            rBox.lengths(),
            std::span<Scalar,cornerCount*Dimension> (boundingVolumes.back().data(), cornerCount * Dimension)
        );
        return true;
    });

    // Compute residuals.
    DynamicArray<Scalar> residual = rhs;
    {
        Eigen::Map<const Eigen::SparseMatrix<Scalar,Eigen::RowMajor,int>> lhsAdaptor(
            lhs.rowCount,
            lhs.columnCount,
            lhs.entries.size(),
            lhs.rowExtents.data(),
            lhs.columnIndices.data(),
            lhs.entries.data());
        Eigen::Map<const Eigen::Matrix<Scalar,Eigen::Dynamic,1>> solutionAdaptor(
            solution.data(),
            solution.size());
        Eigen::Map<Eigen::Matrix<Scalar,Eigen::Dynamic,1>> residualAdaptor(
            residual.data(),
            residual.size());
        residualAdaptor.noalias() -= lhsAdaptor * solutionAdaptor;
    }

    // Collect state samples.
    DynamicArray<Sample> samples(postprocessResolution * postprocessResolution);

    {
        auto logBlock = utils::LoggerSingleton::get().newBlock("scatter postprocess");
        constexpr Scalar epsilon = 1e-10;
        const Scalar postprocessDelta  = (1.0 - 2 * epsilon) / (postprocessResolution - 1);

        mp::ParallelFor<std::size_t>(rThreads).operator()(
            intPow(postprocessResolution, 2),
            [&samples, &solution, &rhs, &residual, &rAssembler, &rMesh, &rBVH, &contiguousCellData, postprocessResolution, postprocessDelta](
                    const std::size_t iSample) -> void {
                const std::size_t iSampleY = iSample / postprocessResolution;
                const std::size_t iSampleX = iSample % postprocessResolution;
                auto& rSample = samples[iSample];
                rSample.position = {epsilon + iSampleX * postprocessDelta,
                                    epsilon + iSampleY * postprocessDelta};

                // Find which cell the global point lies in.
                const auto iMaybeCellData = rBVH.makeView().find(
                    Kernel<Dimension,Scalar>::decayView(rSample.position),
                    std::span<const CellData>(contiguousCellData));

                if (iMaybeCellData != contiguousCellData.size()) {
                    Ref<const CellData> rCellData = contiguousCellData[iMaybeCellData];
                    rSample.cellID = rCellData.id();

                    // Compute sample point in the cell's local space.
                    StaticArray<CellData::ParametricCoordinate,Dimension> localSamplePoint;
                    rCellData.transform(
                        Kernel<Dimension,Scalar>::view(rSample.position),
                        Kernel<Dimension,Scalar>::view(localSamplePoint));

                    // Evaluate the cell's ansatz functions at the local sample point.
                    auto ansatzSpace = rMesh.data().ansatzSpace();
                    std::array<Scalar,Ansatz::size()> ansatzBuffer;
                    ansatzSpace.evaluate(Kernel<Dimension,Scalar>::decayView(
                        localSamplePoint),
                        ansatzBuffer);

                    // Find the entries of the cell's DoFs in the global state vector.
                    const auto& rGlobalIndices = rAssembler[rCellData.id()];

                    // Compute state as an indirect inner product of the solution vector
                    // and the ansatz function values at the local corner coordinates.
                    rSample.state = 0.0;
                    rSample.load = 0.0;
                    rSample.residual = 0.0;
                    for (unsigned iFunction=0u; iFunction<rGlobalIndices.size(); ++iFunction) {
                        rSample.state += solution[rGlobalIndices[iFunction]] * ansatzBuffer[iFunction];
                        rSample.load += rhs[rGlobalIndices[iFunction]] * ansatzBuffer[iFunction];
                        rSample.residual += residual[rGlobalIndices[iFunction]] * ansatzBuffer[iFunction];
                    } // for iFunction in range(rGlobalIndices.size())
                } else {
                    // Could not find the cell that contains the current sample point.
                    rSample.state = NAN;
                }
            });
    }

    // Write XDMF.
    #ifdef CIE_ENABLE_HDF5
        #define CIE_HEAVY_DATA_FORMAT R"("HDF")"
    #else
        #define CIE_HEAVY_DATA_FORMAT R"("XML")"
    #endif

    output << R"(
<Xdmf Version="3.0">
    <Domain>
        <Grid name="root" GridType="Collection" CollectionType="Spatial">

            <Grid Name="boundary" GridType="Uniform">
                <Topology TopologyType="Polyline" NumberOfElements=")" << boundarySegments.size() << R"(" NodesPerElement="2">
                    <DataItem Format=)" << CIE_HEAVY_DATA_FORMAT << R"( Dimensions=")" << boundarySegments.size() << R"( 2">
                        )";
    {
        DynamicArray<std::size_t> data;
        data.reserve(2 * boundarySegments.size());
        std::ranges::for_each(
            std::views::iota(0ul, 2 * boundarySegments.size()),
            [&data] (std::size_t i) {data.push_back(i);});
        output.write(
            std::span<const decltype(data)::value_type>(
                data.data(),
                data.size()),
            "boundaryTopology",
            {static_cast<int>(boundarySegments.size()), 2});
    }
    output << R"(
                    </DataItem>
                </Topology>

                <Geometry Type="XYZ">
                    <DataItem ItemType="Uniform" Format=)" << CIE_HEAVY_DATA_FORMAT << R"( Dimensions=")" << 2 * boundarySegments.size() << R"( 3">
                        )";
    {
        DynamicArray<Scalar> data;
        data.reserve(6 * boundarySegments.size());
        for (const auto& rSegment : boundarySegments) {
            const std::array<Scalar,6> augmentedSegment {
                rSegment[0],
                rSegment[1],
                0.0,
                rSegment[2],
                rSegment[3],
                0.0
            };
            data.insert(
                data.end(),
                augmentedSegment.begin(),
                augmentedSegment.end());
        }
        output.write(
            std::span<const decltype(data)::value_type>(
                data.data(),
                data.size()),
            "boundaryGeometry",
            {static_cast<int>(2 * boundarySegments.size()), 3});
    }
    output << R"(
                    </DataItem>
                </Geometry>
            </Grid>

            <Grid Name="bvh" GridType="Uniform">
                <Topology TopologyType="Quadrilateral" NumberOfElements=")" << boundingVolumes.size() << R"(" NodesPerElement="4">
                    <DataItem Format=)" << CIE_HEAVY_DATA_FORMAT << R"( Dimensions=")" << boundingVolumes.size() << R"( 4">
                        )";
    {
        DynamicArray<int> data;
        data.reserve(boundingVolumes.size() * 4);
        std::ranges::for_each(
            std::views::iota(0ul, boundingVolumes.size()),
            [&data] (int iBoundingBox) {
                const int iPointBegin = 4 * iBoundingBox;
                const std::array<int,4> buffer {
                    iPointBegin,
                    iPointBegin + 1,
                    iPointBegin + 3,
                    iPointBegin + 2
                };
                data.insert(
                    data.end(),
                    buffer.begin(),
                    buffer.end());
            });
        output.write(
            std::span<const decltype(data)::value_type>(
                data.data(),
                data.size()),
            "bvhTopology",
            {static_cast<int>(boundingVolumes.size()), 4});
    }
    output << R"(
                    </DataItem>
                </Topology>

                <Geometry Type="XYZ">
                    <DataItem ItemType="Uniform" Format=)" << CIE_HEAVY_DATA_FORMAT << R"( Dimensions=")" << boundingVolumes.size() * 4 << R"( 3">
                        )";
    {
        DynamicArray<Scalar> data;
        data.reserve(4 * 3 * boundingVolumes.size());
        for (const auto& rBoundingVolume : boundingVolumes) {
            data.push_back(rBoundingVolume[0]);
            data.push_back(rBoundingVolume[1]);
            data.push_back(0.0);
            data.push_back(rBoundingVolume[2]);
            data.push_back(rBoundingVolume[3]);
            data.push_back(0.0);
            data.push_back(rBoundingVolume[4]);
            data.push_back(rBoundingVolume[5]);
            data.push_back(0.0);
            data.push_back(rBoundingVolume[6]);
            data.push_back(rBoundingVolume[7]);
            data.push_back(0.0);
        }
        output.write(
            std::span<const decltype(data)::value_type>(
                data.data(),
                data.size()),
            "bvhGeometry",
            {static_cast<int>(4 * boundingVolumes.size()), 3});
    }
    output << R"(
                    </DataItem>
                </Geometry>
            </Grid>

            <Grid Name="mesh" GridType="Uniform">
                <Topology TopologyType="Quadrilateral" NumberOfElements=")" << rMesh.vertices().size() << R"(" NodesPerElement="4">
                    <DataItem Format=)" << CIE_HEAVY_DATA_FORMAT << R"( Dimensions=")" << rMesh.vertices().size() << R"( 4">
                        )";
    {
        DynamicArray<int> data;
        data.reserve(rMesh.vertices().size() * 4);
        std::ranges::for_each(
            std::views::iota(0ul, rMesh.vertices().size()),
            [&data] (std::size_t iCell) {
                const int i = iCell * 4;
                const std::array<int,4> buffer {
                    i,
                    i + 1,
                    i + 3,
                    i + 2};
                data.insert(
                    data.end(),
                    buffer.begin(),
                    buffer.end());
            });
        output.write(
            std::span<const decltype(data)::value_type>(
                data.data(),
                data.size()),
            "meshTopology",
            {static_cast<int>(rMesh.vertices().size()), 4});
    }
    output << R"(
                    </DataItem>
                </Topology>

                <Geometry Type="XYZ">
                    <DataItem Format=)" << CIE_HEAVY_DATA_FORMAT << R"( Dimensions=")" << 4 * rMesh.vertices().size() << R"( 3">
                        )";
    {
        DynamicArray<Scalar> data;
        const std::array<std::array<CellData::ParametricCoordinate,Dimension>,intPow(2,Dimension)> localCorners {
            std::array<CellData::ParametricCoordinate,Dimension> {-1.0, -1.0},
            std::array<CellData::ParametricCoordinate,Dimension> { 1.0, -1.0},
            std::array<CellData::ParametricCoordinate,Dimension> {-1.0,  1.0},
            std::array<CellData::ParametricCoordinate,Dimension> { 1.0,  1.0}};
        data.reserve(3 * localCorners.size() * rMesh.vertices().size());

        for (const auto& rCell : rMesh.vertices()) {
            for (const auto& rLocalPoint : localCorners) {
                StaticArray<CellData::PhysicalCoordinate,Dimension> globalPoint;
                rCell.data().transform(
                    Kernel<Dimension,Scalar>::view(rLocalPoint),
                    Kernel<Dimension,Scalar>::view(globalPoint));
                data.insert(
                    data.end(),
                    globalPoint.begin(),
                    globalPoint.end());
                data.push_back(0.0);
            } // for rLocalPoint : localCorners
        } // for rCell in rMesh.vertices()

        output.write(
            std::span<const decltype(data)::value_type>(
                data.data(),
                data.size()),
            "meshGeometry",
            {static_cast<int>(4 * rMesh.vertices().size()), 3});
    }
    output << R"(
                    </DataItem>
                </Geometry>

                <Attribute Name="state" Center="Node" AttributeType="Scalar">
                    <DataItem Format=)" << CIE_HEAVY_DATA_FORMAT << R"( Dimensions=")" << 4 * rMesh.vertices().size() << R"(">
                        )";
    {
        DynamicArray<Scalar> data;
        const std::array<std::array<Scalar,Dimension>,intPow(2,Dimension)> localCorners {
            std::array<Scalar,Dimension> {-1.0, -1.0},
            std::array<Scalar,Dimension> { 1.0, -1.0},
            std::array<Scalar,Dimension> {-1.0,  1.0},
            std::array<Scalar,Dimension> { 1.0,  1.0}};

        data.reserve(4 * rMesh.vertices().size());
        DynamicArray<Scalar> ansatzBuffer(rMesh.data().ansatzSpace().size());

        for (const auto& rCell : rMesh.vertices()) {
            const auto& rGlobalIndices = rAssembler[rCell.id()];
            const auto& rAnsatzSpace = rMesh.data().ansatzSpace();
            ansatzBuffer.resize(rAnsatzSpace.size());

            for (const auto& rLocalPoint : localCorners) {
                // Compute ansatz function values at the current local corner.
                rAnsatzSpace.evaluate(rLocalPoint, ansatzBuffer);

                Scalar state = 0.0;
                for (std::size_t iFunction=0; iFunction<ansatzBuffer.size(); ++iFunction) {
                    state += solution[rGlobalIndices[iFunction]] * ansatzBuffer[iFunction];
                }

                data.push_back(state);
            } // for rLocalPoint : localCorners
        } // for rCell in rMesh.vertices()

        output.write(
            std::span<const decltype(data)::value_type>(
                data.data(),
                data.size()),
            "meshState",
            {static_cast<int>(4 * rMesh.vertices().size())});
    }
    output << R"(
                    </DataItem>
                </Attribute>
                <Attribute Name="DoFIDs" Center="Cell" AttributeType="Matrix">
                    <DataItem Format=)" << CIE_HEAVY_DATA_FORMAT << R"( Dimensions=")" << rMesh.vertices().size() << " " << rMesh.data().ansatzSpace().size() << R"(">
                        )";
    {
        DynamicArray<std::size_t> data;

        data.reserve(rMesh.data().ansatzSpace().size() * rMesh.vertices().size());
        DynamicArray<std::size_t> ansatzBuffer(rMesh.data().ansatzSpace().size());

        for (const auto& rCell : rMesh.vertices()) {
            const auto& rGlobalIndices = rAssembler[rCell.id()];
            std::copy(
                rGlobalIndices.begin(),
                rGlobalIndices.end(),
                std::back_inserter(data));
        } // for rCell in rMesh.vertices()

        output.write(
            std::span<const decltype(data)::value_type>(
                data.data(),
                data.size()),
            "DoFIDs",
            {static_cast<int>(rMesh.vertices().size()), static_cast<int>(rMesh.data().ansatzSpace().size())});
    }
    output << R"(
                    </DataItem>
                </Attribute>
                <Attribute Name="axes" Center="Cell" AttributeType="Matrix">
                    <DataItem Format=)" << CIE_HEAVY_DATA_FORMAT << R"( Dimensions=")" << rMesh.vertices().size() << " " << Dimension << R"(">
                        )";
    {
        DynamicArray<std::int8_t> data;

        data.reserve(Dimension * rMesh.vertices().size());

        for (const auto& rCell : rMesh.vertices()) {
            const auto axes = rCell.data().axes();
            for (std::uint8_t iDimension=0; iDimension<Dimension; ++iDimension) {
                const BoundaryID value = axes[iDimension];
                data.push_back(value.getDirection() ? (value.getDimension() + 1) : -(value.getDimension() + 1));
            }
        } // for rCell in rMesh.vertices()

        output.write(
            std::span<const decltype(data)::value_type>(
                data.data(),
                data.size()),
            "cellAxes",
            {static_cast<int>(rMesh.vertices().size()), static_cast<int>(Dimension)});
    }
    output << R"(
                    </DataItem>
                </Attribute>
                <Attribute Name="load" Center="Node" AttributeType="Scalar">
                    <DataItem Format=)" << CIE_HEAVY_DATA_FORMAT << R"( Dimensions=")" << 4 * rMesh.vertices().size() << R"(">
                        )";
    {
        DynamicArray<Scalar> data;
        const std::array<std::array<Scalar,Dimension>,intPow(2,Dimension)> localCorners {
            std::array<Scalar,Dimension> {-1.0, -1.0},
            std::array<Scalar,Dimension> { 1.0, -1.0},
            std::array<Scalar,Dimension> {-1.0,  1.0},
            std::array<Scalar,Dimension> { 1.0,  1.0}};

        data.reserve(4 * rMesh.vertices().size());
        DynamicArray<Scalar> ansatzBuffer(rMesh.data().ansatzSpace().size());

        for (const auto& rCell : rMesh.vertices()) {
            const auto& rGlobalIndices = rAssembler[rCell.id()];
            const auto& rAnsatzSpace = rMesh.data().ansatzSpace();
            ansatzBuffer.resize(rAnsatzSpace.size());

            for (const auto& rLocalPoint : localCorners) {
                // Compute ansatz function values at the current local corner.
                rAnsatzSpace.evaluate(rLocalPoint, ansatzBuffer);

                Scalar state = 0.0;
                for (std::size_t iFunction=0; iFunction<ansatzBuffer.size(); ++iFunction) {
                    state += rhs[rGlobalIndices[iFunction]] * ansatzBuffer[iFunction];
                }

                data.push_back(state);
            } // for rLocalPoint : localCorners
        } // for rCell in rMesh.vertices()

        output.write(
            std::span<const decltype(data)::value_type>(
                data.data(),
                data.size()),
            "meshLoad",
            {static_cast<int>(4 * rMesh.vertices().size())});
    }
    output << R"(
                    </DataItem>
                </Attribute>
                <Attribute Name="residual" Center="Node" AttributeType="Scalar">
                    <DataItem Format=)" << CIE_HEAVY_DATA_FORMAT << R"( Dimensions=")" << 4 * rMesh.vertices().size() << R"(">
                        )";
    {
        DynamicArray<Scalar> data;
        const std::array<std::array<Scalar,Dimension>,intPow(2,Dimension)> localCorners {
            std::array<Scalar,Dimension> {-1.0, -1.0},
            std::array<Scalar,Dimension> { 1.0, -1.0},
            std::array<Scalar,Dimension> {-1.0,  1.0},
            std::array<Scalar,Dimension> { 1.0,  1.0}};

        data.reserve(4 * rMesh.vertices().size());
        DynamicArray<Scalar> ansatzBuffer(rMesh.data().ansatzSpace().size());

        for (const auto& rCell : rMesh.vertices()) {
            const auto& rGlobalIndices = rAssembler[rCell.id()];
            const auto& rAnsatzSpace = rMesh.data().ansatzSpace();
            ansatzBuffer.resize(rAnsatzSpace.size());

            for (const auto& rLocalPoint : localCorners) {
                // Compute ansatz function values at the current local corner.
                rAnsatzSpace.evaluate(rLocalPoint, ansatzBuffer);

                Scalar state = 0.0;
                for (std::size_t iFunction=0; iFunction<ansatzBuffer.size(); ++iFunction) {
                    state += residual[rGlobalIndices[iFunction]] * ansatzBuffer[iFunction];
                }

                data.push_back(state);
            } // for rLocalPoint : localCorners
        } // for rCell in rMesh.vertices()

        output.write(
            std::span<const decltype(data)::value_type>(
                data.data(),
                data.size()),
            "meshResidual",
            {static_cast<int>(4 * rMesh.vertices().size())});
    }
    output << R"(
                    </DataItem>
                </Attribute>
            </Grid>

            <Grid Name="sample" GridType="Uniform">
                <Topology TopologyType="Quadrilateral" NumberOfElements=")" << intPow(postprocessResolution - 1, 2) << R"(" NodePerElement="4">
                    <DataItem Format=)" << CIE_HEAVY_DATA_FORMAT << R"( Dimensions=")" << intPow(postprocessResolution - 1, 2) << R"( 4">
                        )";
    {
        DynamicArray<int> data;
        data.reserve(4 * intPow(postprocessResolution - 1, 2));
        for (unsigned iSampleCellY=0u; iSampleCellY<postprocessResolution - 1; ++iSampleCellY) {
            for (unsigned iSampleCellX=0u; iSampleCellX<postprocessResolution - 1; ++iSampleCellX) {
                const unsigned iBase = iSampleCellY * postprocessResolution + iSampleCellX;
                const std::array<int,4> buffer {
                    static_cast<int>(iBase),
                    static_cast<int>(iBase) + 1,
                    static_cast<int>(iBase + postprocessResolution) + 1,
                    static_cast<int>(iBase + postprocessResolution)};
                data.insert(
                    data.end(),
                    buffer.begin(),
                    buffer.end());
            } // for iSampleCellX in range(proceprocessResolution)
        } // for iSampleCellY in range(proceprocessResolution)
        output.write(
            std::span<const decltype(data)::value_type>(
                data.data(),
                data.size()),
            "sampleTopology",
            {static_cast<int>(intPow(postprocessResolution - 1, 2)), 4});
    }
    output << R"(
                    </DataItem>
                </Topology>

                <Geometry Type="XYZ">
                    <DataItem Format=)" << CIE_HEAVY_DATA_FORMAT << R"( Dimensions=")" << samples.size() << R"( 3">
                        )";
    {
        DynamicArray<Scalar> data;
        data.reserve(3 * samples.size());
        for (const auto& rSample : samples) {
            data.push_back(rSample.position.front());
            data.push_back(rSample.position.back());
            data.push_back(0.0);
        }
        output.write(
            std::span<const decltype(data)::value_type>(
                data.data(),
                data.size()),
            "sampleGeometry",
            {static_cast<int>(samples.size()), 3});
    }
    output << R"(
                    </DataItem>
                </Geometry>

                <Attribute Name="state" Center="Node" AttributeType="Scalar">
                    <DataItem Format=)" << CIE_HEAVY_DATA_FORMAT << R"( Dimensions=")" << samples.size() << R"( 1">
                        )";
    {
        DynamicArray<Scalar> data;
        data.reserve(samples.size());
        std::transform(
            samples.begin(),
            samples.end(),
            std::back_inserter(data),
            [] (const auto& rSample) {
                return rSample.state;
            });
        output.write(
            std::span<const decltype(data)::value_type>(
                data.data(),
                data.size()),
            "sampleState",
            {static_cast<int>(samples.size()), 1});
    }
    output << R"(
                    </DataItem>
                </Attribute>
                <Attribute Name="load" Center="Node" AttributeType="Scalar">
                    <DataItem Format=)" << CIE_HEAVY_DATA_FORMAT << R"( Dimensions=")" << samples.size() << R"( 1">
                        )";
    {
        DynamicArray<Scalar> data;
        data.reserve(samples.size());
        std::transform(
            samples.begin(),
            samples.end(),
            std::back_inserter(data),
            [] (const auto& rSample) {
                return rSample.load;
            });
        output.write(
            std::span<const decltype(data)::value_type>(
                data.data(),
                data.size()),
            "sampleLoad",
            {static_cast<int>(samples.size()), 1});
    }
    output << R"(
                    </DataItem>
                </Attribute>
                <Attribute Name="residual" Center="Node" AttributeType="Scalar">
                    <DataItem Format=)" << CIE_HEAVY_DATA_FORMAT << R"( Dimensions=")" << samples.size() << R"( 1">
                        )";
    {
        DynamicArray<Scalar> data;
        data.reserve(samples.size());
        std::transform(
            samples.begin(),
            samples.end(),
            std::back_inserter(data),
            [] (const auto& rSample) {
                return rSample.residual;
            });
        output.write(
            std::span<const decltype(data)::value_type>(
                data.data(),
                data.size()),
            "sampleResidual",
            {static_cast<int>(samples.size()), 1});
    }
    output << R"(
                    </DataItem>
                </Attribute>

                <Attribute Name="cellID" Center="Node" AttributeType="Scalar">
                    <DataItem Format=)" << CIE_HEAVY_DATA_FORMAT << R"( Dimensions=")" << samples.size() << R"( 1">
                        )";
    {
        DynamicArray<int> data;
        data.reserve(samples.size());
        std::transform(
            samples.begin(),
            samples.end(),
            std::back_inserter(data),
            [] (const auto& rSample) {
                return rSample.cellID;
            });
        output.write(
            std::span<const decltype(data)::value_type>(
                data.data(),
                data.size()),
            "sampleCellID",
            {static_cast<int>(samples.size()), 1});
    }
    output << R"(
                    </DataItem>
                </Attribute>
            </Grid>
        </Grid>
    </Domain>
</Xdmf>
)";
    CIE_END_EXCEPTION_TRACING
}


} // namespace cie::fem
