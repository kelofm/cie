// --- Internal Includes ---
#include "embeddedPoisson2D/definitions.hpp"
#include "embeddedPoisson2D/configuration.hpp"
#include "embeddedPoisson2D/MeshData.hpp"
#include "embeddedPoisson2D/CellData.hpp"
#include "embeddedPoisson2D/mesh.hpp"
#include "embeddedPoisson2D/integration.hpp"
#include "embeddedPoisson2D/constraints.hpp"
#include "embeddedPoisson2D/solver.hpp"
#include "embeddedPoisson2D/postprocessing.hpp"

// --- FEM Includes ---
#include "packages/graph/inc/Assembler.hpp"

// --- GEO Includes ---
#include "packages/io/inc/STLIO.hpp"

// --- Utility Includes ---
#include "packages/commandline/inc/ArgParse.hpp"
#include "packages/logging/inc/LoggerSingleton.hpp"
#include "packages/logging/inc/LogBlock.hpp"
#include "packages/io/inc/json.hpp"


namespace cie::fem {


void readTesselatedDomain(
    Ref<std::istream> rStream,
    Ref<std::vector<Scalar>> rOutput,
    std::span<Scalar,Dimension> bboxBase,
    std::span<Scalar,Dimension> bboxLengths) {
        CIE_BEGIN_EXCEPTION_TRACING
            std::array<Scalar,Dimension> bboxEnd;
            std::transform(
                bboxBase.begin(),
                bboxBase.end(),
                bboxLengths.begin(),
                bboxEnd.begin(),
                std::plus<Scalar>());

            // Read the provided STL input.
            cie::io::STLIO::Input<Scalar,Dimension> io(rStream);
            rOutput.resize(io.triangleCount() * 3 * Dimension);
            io.execute(rOutput);

            // Extend the bounding box.
            std::cout << std::endl;
            for (std::size_t iVertex=0ul; iVertex<3*io.triangleCount(); ++iVertex) {
                for (unsigned iDimension=0u; iDimension<Dimension; ++iDimension) {
                    bboxBase[iDimension] = std::min<Scalar>(
                        bboxBase[iDimension],
                        rOutput[iVertex * Dimension + iDimension]);
                    bboxEnd[iDimension] = std::max<Scalar>(
                        bboxEnd[iDimension],
                        rOutput[iVertex * Dimension + iDimension]);
                }
            }

            // Set the bounding box.
            std::transform(
                bboxEnd.begin(),
                bboxEnd.end(),
                bboxBase.begin(),
                bboxLengths.begin(),
                std::minus<Scalar>());
        CIE_END_EXCEPTION_TRACING
}


int main(Ref<const cie::io::JSONObject> rConfiguration) {
    Mesh mesh;
    mp::ThreadPoolBase threads;

    // Read the dirichlet input and set mesh boundaries.
    const auto tesselatedBoundary = makeBoundary(rConfiguration["dirichlet-1d"]);
    std::array<Scalar,2> bboxBase {
            std::numeric_limits<Scalar>::max(),
            std::numeric_limits<Scalar>::max()},
        bboxLengths{
            std::numeric_limits<Scalar>::lowest(),
            std::numeric_limits<Scalar>::lowest()};

    {
        std::array<Scalar,2> bboxEnd = bboxLengths;
        for (const auto& rSegment : tesselatedBoundary) {
            for (unsigned iPoint=0u; iPoint<2; ++iPoint) {
                for (unsigned iDimension=0u; iDimension<Dimension; ++iDimension) {
                    bboxBase[iDimension] = std::min<Scalar>(
                        bboxBase[iDimension],
                        rSegment[iDimension + iPoint * 2]);
                    bboxEnd[iDimension] = std::max<Scalar>(
                        bboxEnd[iDimension],
                        rSegment[iDimension + iPoint * 2]);
                } // for iDimension
            } // for iPoint
        } // for rSegment
        std::transform(
            bboxEnd.begin(),
            bboxEnd.end(),
            bboxBase.begin(),
            bboxLengths.begin(),
            std::minus<Scalar>());
    }

    // Read domain input and extend mesh boundaries if necessary.
    std::vector<std::pair<
        MeshData::DomainData,
        std::vector<Scalar>
    >> domainTriangles;
    std::vector<std::pair<MeshData::DomainData,Scalar>> domainMap;
    domainMap.emplace_back(
        domainMap.size(),
        0);

    for (const auto& domainConfiguration : rConfiguration["domains"]) {
        const std::string domainType = domainConfiguration["type"].as<std::string>();
        const MeshData::DomainData domainID = domainMap.size();

        // Current settings define the default (fictitious) domain.
        if (domainType == "default") {
            domainMap[0].second = domainConfiguration["scale"].as<double>();
        }

        // Current settings define a triangulated domain.
        else if (domainType == "tesselated") {
            const std::filesystem::path domainFilePath = domainConfiguration["file-path"].as<std::string>();
            auto logBlock = cie::utils::LoggerSingleton::get().newBlock(std::format(
                "reading {}",
                domainFilePath.string()));
            domainTriangles.push_back({
                domainID,
                {}});
            std::ifstream stlFile(
                domainFilePath,
                std::ios::binary);
            readTesselatedDomain(
                stlFile,
                domainTriangles.back().second,
                bboxBase,
                bboxLengths);
            domainMap.emplace_back(
                domainID,
                domainConfiguration["scale"].as<double>());
        } else CIE_THROW(Exception, std::format(
            "unsupported domain type: '{}'",
            domainType))
    }

    // Fill the mesh with cells and boundaries.
    generateMesh(
        mesh,
        bboxBase,
        bboxLengths,
        std::move(domainTriangles),
        domainMap,
        rConfiguration["discretization"]);

    // Find ansatz functions that coincide on opposite boundaries.
    // In adjacent cells, these ansatz functions will have to map
    // to the same DoF in the assembled system.
    const auto ansatzMap = makeAnsatzMap(
        mesh.data().ansatz(0ul),
        /*integrationOrder=*/5,
        utils::Comparison<Scalar>(
            /*absoluteTolerance =*/ std::numeric_limits<Scalar>::epsilon(),
            /*relativeTolerance =*/ 1e2 * std::numeric_limits<Scalar>::epsilon()));

    // Build a factory detects the topology of the mesh and issues indices to DoFs.
    Assembler assembler;

    {
        auto logBlock = utils::LoggerSingleton::get().newBlock("parse mesh topology");
        assembler.addGraph<Cell>(
            mesh,
            [&mesh] (VertexID id) -> Ref<const Cell> {
                const auto& rCells = mesh.data().cells();
                return *std::lower_bound(
                    rCells.begin(),
                    rCells.end(),
                    id,
                    [] (Ref<const Cell> rCell, VertexID id) {
                        return rCell.id() < id;});
            },
            ansatzMap,
            mesh.data().ansatz(0ul).size());
    } // parse mesh topology

    // Construct a bounding volume hierarchy over the cells to accelerate
    // point membership tests. Running on accelerator devices also requires
    // - the cell data to be available in a contiguous array
    // - the cell data to be self contained (no pointers and heap storage)
    auto bvh = makeBoundingVolumeHierarchy(
        mesh.data().cells(),
        bboxBase,
        bboxLengths);
    const auto bvhView = bvh.makeView();

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
            entries,
            threads);
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
        mesh.data().cells(),
        assembler,
        lhs,
        rConfiguration["discretization"]["integration"],
        threads);

    const auto boundarySegments = imposeBoundaryConditions(
        mesh,
        tesselatedBoundary,
        assembler,
        bvhView,
        mesh.data().cells(),
        lhs,
        rhs,
        rConfiguration["dirichlet-1d"]);

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
            rConfiguration["linear-system"]);
    } // solve

    {
        auto logBlock = utils::LoggerSingleton::get().newBlock("postprocess");
        postprocess(
            bboxBase,
            bboxLengths,
            lhs,
            solution,
            rhs,
            mesh,
            bvh,
            assembler,
            rConfiguration["discretization"]["postprocessing"],
            threads);
    }

    return 0;
}


} // namespace cie::fem

int main(int argc, const char** argv) {
    const auto logBlock = cie::utils::LoggerSingleton::get().newBlock("main");
    cie::utils::ArgParse parser("2D Poisson Example");
    parser
        .addPositional("config-path")
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

    const std::filesystem::path configPath = arguments.get<std::filesystem::path>("config-path");
    cie::io::JSONObject configuration;

    // Read the provided configuration file.
    CIE_BEGIN_EXCEPTION_TRACING
        std::cout << "Read input configuration from " << configPath << ".\n";
        std::ifstream configFile(configPath);
        configuration = cie::io::JSONObject(configFile);
    CIE_END_EXCEPTION_TRACING

    // Validate the provided configuration and apply defaults to it.
    CIE_BEGIN_EXCEPTION_TRACING
        cie::io::JSONSchema configSchema;
        cie::fem::makeSchema(configSchema);
        configSchema.validateAndFillDefaults(configuration);
    CIE_END_EXCEPTION_TRACING

    // Write the augmented configuration.
    CIE_BEGIN_EXCEPTION_TRACING
        const std::filesystem::path appliedConfigPath =
            configPath.parent_path() / std::format(
                "{}_applied{}",
                std::filesystem::path::string_type(configPath.stem()),
                std::filesystem::path::string_type(configPath.extension()));
        std::cout << "Write applied configuration to " << appliedConfigPath << ".\n";
        std::ofstream appliedConfigFile(appliedConfigPath);
        configuration.prettyPrint(appliedConfigFile);
    CIE_END_EXCEPTION_TRACING

    try {
        return cie::fem::main(configuration);
    } catch (cie::Exception& rException) {
        std::cerr << rException.what() << std::endl;
        return 1;
    }
}
