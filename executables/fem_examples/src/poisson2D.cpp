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
#include "poisson2D/xdmf.hpp"

// --- FEM Includes ---
#include "packages/graph/inc/connectivity.hpp"
#include "packages/graph/inc/Assembler.hpp"


namespace cie::fem {


int main(Ref<const utils::ArgParse::Results> rArguments) {
    Mesh mesh;
    mp::ThreadPoolBase threads;

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

    // Compute element contributions and assemble them into the matrix
    {
        integrateStiffness(
            mesh,
            assembler,
            lhs,
            rArguments,
            threads);
    }

    const auto boundarySegments = imposeBoundaryConditions(
        mesh,
        assembler,
        bvhView,
        contiguousCellData,
        lhs,
        rhs,
        rArguments);

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

    postprocess(
        lhs,
        solution,
        rhs,
        mesh,
        contiguousCellData,
        boundarySegments,
        bvh,
        assembler,
        rArguments,
        threads);

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
