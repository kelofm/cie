// --- External Includes ---
#include "Eigen/Sparse"
#include "Eigen/IterativeLinearSolvers"
#include "Eigen/src/SparseCholesky/SimplicialCholesky.h"

// --- Internal Includes ---
#include "packages/commandline/inc/ArgParse.hpp"
#include "poisson2D/definitions.hpp"
#include "poisson2D/MeshData.hpp"
#include "poisson2D/CellData.hpp"
#include "poisson2D/BoundaryData.hpp"
#include "poisson2D/mesh.hpp"
#include "poisson2D/integration.hpp"
#include "poisson2D/constraints.hpp"

// --- FEM Includes ---
#include "packages/graph/inc/connectivity.hpp"
#include "packages/graph/inc/Assembler.hpp"

// --- STL Includes ---
#include <ranges> // ranges::iota


namespace cie::fem {



int main(Ref<const utils::ArgParse::Results> rArguments) {
    Mesh mesh;
    mp::ThreadPoolBase threads;
    const unsigned postprocessResolution = rArguments.get<std::size_t>("scatter-resolution");

    // Fill the mesh with cells and boundaries.
    generateMesh(mesh, rArguments);

    {
        auto logBlock = utils::LoggerSingleton::get().newBlock("write graphml");
        io::GraphML::GraphML::Output("poisson2D.graphml")(mesh);
    }

    // Find ansatz functions that coincide on opposite boundaries.
    // In adjacent cells, these ansatz functions will have to map
    // to the same DoF in the assembled system.
    const StaticArray<Scalar,5> samples {-1.0, -0.5, 0.0, 0.5, 1.0};
    const auto ansatzMap = makeAnsatzMap(
        mesh.data().ansatzSpace(),
        samples,
        utils::Comparison<Scalar>(
            /*absoluteTolerance =*/ 1e-8,
            /*relativeTolerance =*/ 1e-6));

    // Build a factory detects the topology of the mesh and issues indices to DoFs.
    Assembler assembler;

    {
        auto logBlock = utils::LoggerSingleton::get().newBlock("parse mesh topology");
        assembler.addGraph(
            mesh,
            [&mesh]([[maybe_unused]] Ref<const Mesh::Vertex> rVertex) -> std::size_t {
                const Ansatz& rAnsatz = mesh.data().ansatzSpace();
                return rAnsatz.size();
            },
            [&ansatzMap, &mesh = std::as_const(mesh)](Ref<const Mesh::Edge> rEdge, Assembler::DoFPairIterator it) {
                const auto sourceAxes = mesh.find(rEdge.source()).value().data().axes();
                const auto targetAxes = mesh.find(rEdge.target()).value().data().axes();
                ansatzMap.getPairs(
                    OrientedBoundary<Dimension>(sourceAxes, rEdge.data().boundary),
                    OrientedBoundary<Dimension>(targetAxes, rEdge.data().boundary),
                    it);
            });
    } // parse mesh topology

    // Construct a bounding volume hierarchy over the cells to accelerate
    // point membership tests. Running on accelerator devices also requires
    // - the cell data to be available in a contiguous array
    // - the cell data to be self contained (no pointers and heap storage)
    auto bvh = makeBoundingVolumeHierarchy(mesh);
    const auto bvhView = bvh.makeView();
    DynamicArray<CellData> contiguousCellData(mesh.vertices().size());
    std::ranges::transform(
        mesh.vertices(),
        contiguousCellData.data(),
        [](Ref<const Mesh::Vertex> rCell){return rCell.data();});

    // Create empty CSR matrix
    int rowCount, columnCount;
    DynamicArray<int> rowExtents, columnIndices;
    DynamicArray<double> entries;
    {
        auto logBlock = utils::LoggerSingleton::get().newBlock("compute sparsity pattern");
        assembler.makeCSRMatrix(
            rowCount,
            columnCount,
            rowExtents,
            columnIndices,
            entries);
    }
    DynamicArray<Scalar> rhs(rowCount, 0.0);
    CSRWrapper lhs {
        .rowCount = rowCount,
        .columnCount = columnCount,
        .rowExtents = rowExtents,
        .columnIndices = columnIndices,
        .entries = entries
    };

    const auto boundarySegments = imposeBoundaryConditions(
        mesh,
        assembler,
        bvhView,
        contiguousCellData,
        lhs,
        rhs,
        rArguments);

    // Compute element contributions and assemble them into the matrix
    integrateStiffness(
        mesh,
        assembler,
        lhs,
        rArguments,
        threads);

    // Solve the linear system.
    DynamicArray<Scalar> solution(rhs.size());
    {
        auto logBlock = utils::LoggerSingleton::get().newBlock("solve linear system");
        using EigenSparseMatrix = Eigen::SparseMatrix<Scalar,Eigen::RowMajor,int>;
        Eigen::Map<EigenSparseMatrix> lhsAdaptor(
            rowCount,
            columnCount,
            entries.size(),
            rowExtents.data(),
            columnIndices.data(),
            entries.data());
        Eigen::Map<Eigen::Matrix<Scalar,Eigen::Dynamic,1>> rhsAdaptor(rhs.data(), rhs.size(), 1);
        Eigen::Map<Eigen::Matrix<Scalar,Eigen::Dynamic,1>> solutionAdaptor(solution.data(), solution.size(), 1);

        Eigen::ConjugateGradient<
            EigenSparseMatrix,
            Eigen::Lower | Eigen::Upper,
            Eigen::DiagonalPreconditioner<Scalar>
        > solver;
        solver.setMaxIterations(int(5e3));
        solver.setTolerance(1e-6);
        //Eigen::SimplicialLLT<EigenSparseMatrix> solver;

        solver.compute(lhsAdaptor);
        solutionAdaptor = solver.solve(rhsAdaptor);

        std::cout << solver.iterations() << " iterations "
                  << solver.error()      << " residual\n";
    } // solve

    // XDMF output.
    {
        // Collect output for the bounding volume hierarchy.
        DynamicArray<StaticArray<Scalar,Dimension*intPow(2,Dimension)>> boundingVolumes; // {p0x, p0y, p1x, p1y, p2x, p2y, p3x, p3y}
        bvh.visit([&boundingVolumes](const auto& rBox) -> bool {
            constexpr unsigned cornerCount = intPow(2,Dimension);
            boundingVolumes.emplace_back();
            geo::Box<Dimension,Scalar>::makeCorners(
                rBox.base(),
                rBox.lengths(),
                std::span<Scalar,cornerCount*Dimension> (boundingVolumes.back().data(), cornerCount * Dimension)
            );
            return true;
        });

        // Collect state samples.
        struct Sample {
            StaticArray<CellData::GlobalCoordinate,Dimension>   position;
            Scalar                                              state;
            unsigned                                            cellID;
        }; // struct Sample
        DynamicArray<Sample> samples(postprocessResolution * postprocessResolution);

        {
            auto logBlock = utils::LoggerSingleton::get().newBlock("scatter postprocess");
            constexpr Scalar epsilon = 1e-10;
            const Scalar postprocessDelta  = (1.0 - 2 * epsilon) / (postprocessResolution - 1);

            mp::ParallelFor<unsigned>(threads).operator()(
                intPow(postprocessResolution, 2),
                [&samples, &solution, &assembler, bvhView, &mesh, &contiguousCellData, postprocessResolution, postprocessDelta](
                        const unsigned iSample) -> void {
                    const unsigned iSampleY = iSample / postprocessResolution;
                    const unsigned iSampleX = iSample % postprocessResolution;
                    auto& rSample = samples[iSample];
                    rSample.position = {epsilon + iSampleX * postprocessDelta,
                                        epsilon + iSampleY * postprocessDelta};

                    // Find which cell the global point lies in.
                    const auto iMaybeCellData = bvhView.find(
                        Kernel<Dimension,Scalar>::decayView(rSample.position),
                        std::span<const CellData>(contiguousCellData)
                    );

                    if (iMaybeCellData != contiguousCellData.size()) {
                        Ref<const CellData> rCellData = contiguousCellData[iMaybeCellData];
                        rSample.cellID = rCellData.id();

                        // Compute sample point in the cell's local space.
                        StaticArray<CellData::LocalCoordinate,Dimension> localSamplePoint;
                        rCellData.transform(
                            Kernel<Dimension,Scalar>::view(rSample.position),
                            Kernel<Dimension,Scalar>::view(localSamplePoint));

                        // Evaluate the cell's ansatz functions at the local sample point.
                        auto ansatzSpace = mesh.data().ansatzSpace();
                        std::array<Scalar,Ansatz::size()> ansatzBuffer;
                        ansatzSpace.evaluate(Kernel<Dimension,Scalar>::decayView(
                            localSamplePoint),
                            ansatzBuffer);

                        // Find the entries of the cell's DoFs in the global state vector.
                        const auto& rGlobalIndices = assembler[rCellData.id()];

                        // Compute state as an indirect inner product of the solution vector
                        // and the ansatz function values at the local corner coordinates.
                        rSample.state = 0;
                        for (unsigned iFunction=0u; iFunction<rGlobalIndices.size(); ++iFunction) {
                            rSample.state += solution[rGlobalIndices[iFunction]] * ansatzBuffer[iFunction];
                        } // for iFunction in range(rGlobalIndices.size())
                    } else {
                        // Could not find the cell that contains the current sample point.
                        rSample.state = NAN;
                    }
                });
        }

        // Write XDMF.
        std::ofstream xdmf("poisson2D.xdmf");
        xdmf << R"(
<Xdmf Version="3.0">
    <Domain>
        <Grid name="root" GridType="Collection" CollectionType="Spatial">

            <Grid Name="boundary" GridType="Uniform">
                <Topology TopologyType="Polyline" NumberOfElements=")" << boundarySegments.size() << R"(" NodesPerElement="2">
                    <DataItem Format="XML" Dimensions=")" << boundarySegments.size() << R"( 2">
                        )"; {
                            for (unsigned iSegment=0u; iSegment<boundarySegments.size(); ++iSegment) {
                                xdmf << 2 * iSegment << " " << 2 * iSegment + 1 << "\n                        ";
                            }
                        } xdmf << R"(
                    </DataItem>
                </Topology>

                <Geometry Type="XYZ">
                    <DataItem ItemType="Uniform" Format="XML" Dimensions=")" << 2 * boundarySegments.size() << R"( 3">
                        )"; {
                            for (const auto& rSegment : boundarySegments) {
                                xdmf << rSegment[0] << " " << rSegment[1] << " 0" << "\n                        "
                                     << rSegment[2] << " " << rSegment[3] << " 0" << "\n                        ";
                            }
                        };
                        xdmf << R"(
                    </DataItem>
                </Geometry>

                <Attribute Name="level" Center="Cell" AttributeType="Scalar">
                    <DataItem Format="XML" Dimensions=")" << boundarySegments.size() << R"( 1">
                        )"; {
                            for (const auto& rSegment : boundarySegments) {
                                xdmf << rSegment.back() << "\n                        ";
                            }
                        } xdmf << R"(
                    </DataItem>
                </Attribute>
            </Grid>

            <Grid Name="bvh" GridType="Uniform">
                <Topology TopologyType="Quadrilateral" NumberOfElements=")" << boundingVolumes.size() << R"(" NodesPerElement="4">
                    <DataItem Format="XML" Dimensions=")" << boundingVolumes.size() << R"( 4">
                        )"; {
                            for (unsigned iBoundingVolume=0u; iBoundingVolume<boundingVolumes.size(); ++iBoundingVolume) {
                                const unsigned iPointBegin = 4 * iBoundingVolume;
                                xdmf << iPointBegin << " "
                                     << iPointBegin + 1 << " "
                                     << iPointBegin + 3 << " "
                                     << iPointBegin + 2
                                     << "\n                        ";
                            }
                        } xdmf << R"(
                    </DataItem>
                </Topology>

                <Geometry Type="XYZ">
                    <DataItem ItemType="Uniform" Format="XML" Dimensions=")" << boundingVolumes.size() * 4 << R"( 3">
                        )"; {
                            for (const auto& rBoundingVolume : boundingVolumes) {
                                for (unsigned iCorner=0u; iCorner<4u; ++iCorner) {
                                    for (unsigned iDimension=0u; iDimension<Dimension; ++iDimension)
                                        xdmf << rBoundingVolume[Dimension * iCorner + iDimension] << " ";
                                    xdmf << 0 << "\n                        ";
                                }
                            }
                        } xdmf << R"(
                    </DataItem>
                </Geometry>
            </Grid>

            <Grid Name="mesh" GridType="Uniform">
                <Topology TopologyType="Quadrilateral" NumberOfElements=")" << mesh.vertices().size() << R"(" NodesPerElement="4">
                    <DataItem Format="XML" Dimensions=")" << mesh.vertices().size() << R"( 4">
                        )"; {
                            for (unsigned iCell=0u; iCell<mesh.vertices().size(); ++iCell) {
                                const unsigned i = 4 * iCell;
                                xdmf << i << " " << i + 1 << " " << i + 3 << " " << i + 2 << "\n                        ";
                            } // for iCell in range(mesh.vertices().size())
                        } xdmf << R"(
                    </DataItem>
                </Topology>

                <Geometry Type="XYZ">
                    <DataItem Format="XML" Dimensions=")" << 4 * mesh.vertices().size() << R"( 3">
                        )"; {
                            StaticArray<StaticArray<CellData::LocalCoordinate,Dimension>,intPow(2,Dimension)> localCorners {
                                {-1.0, -1.0},
                                { 1.0, -1.0},
                                {-1.0,  1.0},
                                { 1.0,  1.0}
                            };

                            for (const auto& rCell : mesh.vertices()) {
                                for (const auto& rLocalPoint : localCorners) {
                                    StaticArray<CellData::GlobalCoordinate,Dimension> globalPoint;
                                    rCell.data().transform(
                                        Kernel<Dimension,Scalar>::view(rLocalPoint),
                                        Kernel<Dimension,Scalar>::view(globalPoint));
                                    for (auto c : globalPoint) xdmf << c << " ";
                                    xdmf << "0\n                        ";
                                } // for rLocalPoint : localCorners
                            } // for rCell in mesh.vertices()
                        } xdmf << R"(
                    </DataItem>
                </Geometry>

                <Attribute Name="state" Center="Node" AttributeType="Scalar">
                    <DataItem Format="XML" Dimensions=")" << 4 * mesh.vertices().size() << R"( 1">
                        )"; {
                            StaticArray<StaticArray<Scalar,Dimension>,intPow(2,Dimension)> localCorners {
                                {-1.0, -1.0},
                                { 1.0, -1.0},
                                {-1.0,  1.0},
                                { 1.0,  1.0}
                            };
                            DynamicArray<Scalar> ansatzBuffer(mesh.data().ansatzSpace().size());

                            for (const auto& rCell : mesh.vertices()) {
                                const auto& rGlobalIndices = assembler[rCell.id()];
                                const auto& rAnsatzSpace = mesh.data().ansatzSpace();
                                ansatzBuffer.resize(rAnsatzSpace.size());

                                for (const auto& rLocalPoint : localCorners) {
                                    // Compute ansatz function values at the current local corner.
                                    rAnsatzSpace.evaluate(rLocalPoint, ansatzBuffer);

                                    // Compute state as an indirect inner product of the solution vector
                                    // and the ansatz function values at the local corner coordinates.
                                    Scalar state = 0;
                                    for (unsigned iFunction=0u; iFunction<rGlobalIndices.size(); ++iFunction) {
                                        state += solution[rGlobalIndices[iFunction]] * ansatzBuffer[iFunction];
                                    } // for iFunction in range(rGlobalIndices.size())

                                    xdmf << state << "\n                        ";
                                } // for rLocalPoint : localCorners
                            } // for rCell in mesh.vertices()
                        } xdmf << R"(
                    </DataItem>
                </Attribute>
            </Grid>

            <Grid Name="sample" GridType="Uniform">
                <Topology TopologyType="Quadrilateral" NumberOfElements=")" << intPow(postprocessResolution - 1, 2) << R"(" NodePerElement="4">
                    <DataItem Format="XML" Dimensions=")" << intPow(postprocessResolution - 1, 2) << R"( 4">
                        )"; {
                            for (unsigned iSampleCellY=0u; iSampleCellY<postprocessResolution - 1; ++iSampleCellY) {
                                for (unsigned iSampleCellX=0u; iSampleCellX<postprocessResolution - 1; ++iSampleCellX) {
                                    const unsigned iBase = iSampleCellY * postprocessResolution + iSampleCellX;
                                    xdmf << iBase << " " << iBase + 1 << " " << iBase + postprocessResolution +1 << " " << iBase + postprocessResolution << "\n                        ";
                                } // for iSampleCellX in range(proceprocessResolution)
                            } // for iSampleCellY in range(proceprocessResolution)
                        } xdmf << R"(
                    </DataItem>
                </Topology>

                <Geometry Type="XYZ">
                    <DataItem Format="XML" Dimensions=")" << samples.size() << R"( 3">
                        )"; {
                            for (const auto& rSample : samples) {
                                xdmf << rSample.position[0] << " "
                                     << rSample.position[1]
                                     << " 0\n                        ";
                            } // for rSample in samples
                        } xdmf << R"(
                    </DataItem>
                </Geometry>

                <Attribute Name="state" Center="Node" AttributeType="Scalar">
                    <DataItem Format="XML" Dimensions=")" << samples.size() << R"( 1">
                        )"; {
                            for (const auto& rSample : samples) {
                                xdmf << rSample.state << "\n                        ";
                            } // for rSample in samples
                        } xdmf << R"(
                    </DataItem>
                </Attribute>

                <Attribute Name="cellID" Center="Node" AttributeType="Scalar">
                    <DataItem Format="XML" Dimensions=")" << samples.size() << R"( 1">
                        )"; {
                            for (const auto& rSample : samples) {
                                xdmf << rSample.cellID << "\n                        ";
                            } // for rSample in samples
                        } xdmf << R"(
                    </DataItem>
                </Attribute>
            </Grid>
        </Grid>
    </Domain>
</Xdmf>
)";
    }

    return 0;
}


} // namespace cie::fem

int main(int argc, const char** argv) {
    const auto logBlock = cie::utils::LoggerSingleton::get().newBlock("main");
    cie::utils::ArgParse parser("2D Poisson Example");
    parser
        .addKeyword(
            {"-r", "--resolution"},
            cie::utils::ArgParse::DefaultValue {"31"},
            cie::utils::ArgParse::validatorFactory<std::size_t>(),
            "Resolution (number of nodes per direction.)")
        .addKeyword(
            {"--boundary-resolution"},
            cie::utils::ArgParse::DefaultValue {"20"},
            cie::utils::ArgParse::validatorFactory<std::size_t>(),
            "Number of nodes discretizing the boundary circle.")
        .addKeyword(
            {"--scatter-resolution"},
            cie::utils::ArgParse::DefaultValue {"100"},
            cie::utils::ArgParse::validatorFactory<std::size_t>(),
            "Number of sample points per direction for scatter post-processing.")
        .addKeyword(
            {"--penalty-factor"},
            cie::utils::ArgParse::DefaultValue {"1e-2"},
            cie::utils::ArgParse::validatorFactory<double>(),
            "Penalty value for the weak imposition of Dirichlet boundary conditions.")
        .addKeyword(
            {"--integrand-batch-size"},
            cie::utils::ArgParse::DefaultValue {"0x8000"},
            cie::utils::ArgParse::validatorFactory<std::size_t>(),
            "Number of quadrature points to evaluate at once.")
        .addFlag(
            {"--gpu"},
            "Use the GPU for integration and postprocessing.")
        .addFlag(
            {"-h", "--help"},
            "Print this help and exit.")
        ;

    cie::utils::ArgParse::Results arguments;
    try {
        arguments = parser.parseArguments(argc - 1, argv + 1);
    } catch (cie::Exception& rException) {
        std::cerr << rException.what() << std::endl;
        return 1;
    }

    if (arguments.get<bool>("h")) {
        parser.help(std::cout);
        return 0;
    }

    try {
        return cie::fem::main(arguments);
    } catch (cie::Exception& rException) {
        std::cerr << rException.what() << std::endl;
        return 1;
    }
}
