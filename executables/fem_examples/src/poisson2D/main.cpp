// --- Internal Includes ---
#include "packages/commandline/inc/ArgParse.hpp"
#include "poisson2D/definitions.hpp"
#include "poisson2D/MeshData.hpp"
#include "poisson2D/CellData.hpp"
#include "poisson2D/mesh.hpp"
#include "poisson2D/integration.hpp"
#include "poisson2D/constraints.hpp"
#include "poisson2D/solver.hpp"
#include "poisson2D/postprocessing.hpp"

// --- FEM Includes ---
#include "packages/graph/inc/Assembler.hpp"

// --- Utility Includes ---
#include "packages/logging/inc/LoggerSingleton.hpp"
#include "packages/logging/inc/LogBlock.hpp"


namespace cie::fem {


int main(Ref<const utils::ArgParse::Results> rArguments) {
    Mesh mesh;
    mp::ThreadPoolBase threads;

    // Read the boundary input and set mesh boundaries.
    const auto tesselatedBoundary = makeBoundary(rArguments);
    std::array<Scalar,2> meshBase {
            std::numeric_limits<Scalar>::max(),
            std::numeric_limits<Scalar>::max()},
        meshLengths{
            std::numeric_limits<Scalar>::lowest(),
            std::numeric_limits<Scalar>::lowest()};
    for (const auto& rSegment : tesselatedBoundary) {
        for (unsigned iDimension=0u; iDimension<2; ++iDimension) {
            for (unsigned iPoint=0u; iPoint<2; ++iPoint) {
                meshBase[iDimension] = std::min<Scalar>(
                    meshBase[iDimension],
                    rSegment[iDimension + iPoint * 2]);
                meshLengths[iDimension] = std::max<Scalar>(
                    meshLengths[iDimension],
                    rSegment[iDimension + iPoint * 2]);
            }
        }
    }
    meshLengths.front() = (1.0 + 2e-2) * (meshLengths.front() - meshBase.front());
    meshLengths.back()  = (1.0 + 2e-2) * (meshLengths.back() - meshBase.back());
    meshBase.front() -= 1e-2 / (1.0 + 2e-2) * (meshLengths.front());
    meshBase.back() -= 1e-2 / (1.0 + 2e-2) * (meshLengths.back());
    std::cout << std::format(
        "mesh covers\n\tx in [{}, {}]\n\ty in [{}, {}]\n",
        meshBase.front(), meshBase.front() + meshLengths.front(),
        meshBase.back(), meshBase.back() + meshLengths.back());

    // Fill the mesh with cells and boundaries.
    generateMesh(
        mesh,
        meshBase,
        meshLengths,
        rArguments);

    {
        auto logBlock = utils::LoggerSingleton::get().newBlock("write graphml");
        io::GraphML::GraphML::Output("poisson2D.graphml")(mesh);
    }

    // Find ansatz functions that coincide on opposite boundaries.
    // In adjacent cells, these ansatz functions will have to map
    // to the same DoF in the assembled system.
    const auto ansatzMap = makeAnsatzMap(
        mesh.data().ansatz(0ul),
        /*integrationOrder=*/5,
        utils::Comparison<Scalar>(
            /*absoluteTolerance =*/ 1e-8,
            /*relativeTolerance =*/ 1e-6));

    // Build a factory detects the topology of the mesh and issues indices to DoFs.
    Assembler assembler;

    {
        auto logBlock = utils::LoggerSingleton::get().newBlock("parse mesh topology");
        assembler.addGraph(
            mesh,
            ansatzMap,
            mesh.data().ansatz(0ul).size());
    } // parse mesh topology

    // Construct a bounding volume hierarchy over the cells to accelerate
    // point membership tests. Running on accelerator devices also requires
    // - the cell data to be available in a contiguous array
    // - the cell data to be self contained (no pointers and heap storage)
    auto bvh = makeBoundingVolumeHierarchy(mesh, meshBase, meshLengths);
    const auto bvhView = bvh.makeView();
    DynamicArray<CellData> contiguousCellData(mesh.vertices().size());
    std::ranges::transform(
        mesh.vertices(),
        contiguousCellData.data(),
        [](Ref<const Mesh::Vertex> rCell){return rCell.data();});

    // Create empty CSR matrix
    int rowCount, columnCount;
    DynamicArray<int> rowExtents, columnIndices;
    DynamicArray<Scalar> entries;
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
    linalg::CSRView<Scalar,int> lhs (
        columnCount,
        rowExtents,
        columnIndices,
        entries);

    // Compute element contributions and assemble them into the matrix
    integrateStiffness(
        mesh,
        contiguousCellData,
        assembler,
        lhs,
        rArguments,
        threads);

    const auto boundarySegments = imposeBoundaryConditions(
        mesh,
        tesselatedBoundary,
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
        solve(
            lhs,
            solution,
            rhs,
            assembler,
            threads,
            rArguments);
    } // solve

    {
        auto logBlock = utils::LoggerSingleton::get().newBlock("postprocess");
        postprocess(
            meshBase,
            meshLengths,
            lhs,
            solution,
            rhs,
            mesh,
            contiguousCellData,
            bvh,
            assembler,
            rArguments,
            threads);
    }

    return 0;
}


} // namespace cie::fem

int main(int argc, const char** argv) {
    const auto logBlock = cie::utils::LoggerSingleton::get().newBlock("main");
    cie::utils::ArgParse parser("2D Poisson Example");
    parser
        .addKeyword(
            {"--boundary-file-path"},
            cie::utils::ArgParse::DefaultValue {"boundary.csv"},
            //cie::utils::ArgParse::validatorFactory<std::filesystem::path>(),
            "Path to the boundary definition file.")
        .addKeyword(
            {"-r", "--resolution"},
            cie::utils::ArgParse::DefaultValue {"31"},
            cie::utils::ArgParse::validatorFactory<std::size_t>(),
            "Resolution (number of nodes per direction.)")
        .addKeyword(
            {"--basis"},
            cie::utils::ArgParse::DefaultValue {"legendre"},
            [] (cie::utils::ArgParse::ValueView argument) {
                std::string arg;
                std::copy(argument.begin(), argument.end(), std::back_inserter(arg));
                return arg == "legendre" || arg == "lagrange";
            },
            "Type of basis polynomials to use (lagrange or integrated legendre)")
        .addKeyword(
            {"--min-boundary-tree-depth"},
            cie::utils::ArgParse::DefaultValue {"3"},
            cie::utils::ArgParse::validatorFactory<std::size_t>(),
            "Minimum number of splits before considering boundary segments for integration.")
        .addKeyword(
            {"--max-boundary-tree-depth"},
            cie::utils::ArgParse::DefaultValue {"15"},
            cie::utils::ArgParse::validatorFactory<std::size_t>(),
            "Maximum number of segment splits on the boundary.")
        .addKeyword(
            {"--min-boundary-segment-norm"},
            cie::utils::ArgParse::DefaultValue {"1e-12"},
            cie::utils::ArgParse::validatorFactory<double>(),
            "Minimum size of a boundary segment to integrate.")
        .addKeyword(
            {"--solver"},
            cie::utils::ArgParse::DefaultValue {"cg-eigen"},
            [] (const cie::utils::ArgParse::ValueView& rValue) -> bool {
                const std::string value(rValue.begin(), rValue.end());
                const std::vector<std::string> choices {
                    "cg",
                    "cg-sycl",
                    "cg-eigen",
                    "p-multigrid",
                    "jacobi",
                    "llt-eigen"};
                return std::any_of(
                    choices.begin(),
                    choices.end(),
                    [&value] (cie::Ref<const std::string> rChoice) {return value == rChoice;});
            },
            "Linear solver.")
        .addKeyword(
            {"--reordering"},
            cie::utils::ArgParse::DefaultValue {"none"},
            [] (const cie::utils::ArgParse::ValueView& rValue) -> bool {
                const std::string value(rValue.begin(), rValue.end());
                return value == "none" || value == "cuthill-mckee" || value == "reverse-cuthill-mckee";
            },
            "LHS matrix reordering strategy.")
        .addKeyword(
            {"--scatter-resolution"},
            cie::utils::ArgParse::DefaultValue {"100"},
            cie::utils::ArgParse::validatorFactory<std::size_t>(),
            "Number of sample points per direction for scatter post-processing.")
        .addKeyword(
            {"--penalty-factor"},
            cie::utils::ArgParse::DefaultValue {"1e6"},
            cie::utils::ArgParse::validatorFactory<double>(),
            "Penalty value for the weak imposition of Dirichlet boundary conditions.")
        .addKeyword(
            {"--integrand-batch-size"},
            cie::utils::ArgParse::DefaultValue {"0x8000"},
            cie::utils::ArgParse::validatorFactory<std::size_t>(),
            "Number of quadrature points to evaluate at once.")
        .addFlag(
            {"--sycl"},
            "Use the GPU for integration and postprocessing.")
        .addFlag(
            {"--write-linear-system"},
            "Write the linear system in matrix market format.")
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
