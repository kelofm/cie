// --- Internal Includes ---
#include "poisson2D/definitions.hpp"
#include "poisson2D/mesh.hpp"

// --- FEM Includes ---
#include "packages/maths/inc/LegendrePolynomial.hpp"
#include "packages/maths/inc/LagrangePolynomial.hpp"

// --- Utility Includes ---
#include "packages/logging/inc/LoggerSingleton.hpp"
#include "packages/io/inc/json.hpp"


namespace cie::fem {


void parseGeometricDiscretization(
    Ref<const cie::io::JSONObject> rConfiguration,
    std::span<const Scalar> bboxBase,
    std::span<const Scalar> bboxLengths,
    std::span<Scalar,Dimension> meshBase,
    std::span<Scalar,Dimension> meshLengths,
    std::span<std::size_t,Dimension> meshResolution) {
        CIE_BEGIN_EXCEPTION_TRACING
            const std::string discretizationType = rConfiguration["type"].as<std::string>();

            // The configuration defines a uniform cartesian grid.
            if (discretizationType == "uniform-cartesian-grid") {
                const Scalar meshLength = *std::max_element(bboxLengths.begin(), bboxLengths.end());
                std::fill(
                    meshLengths.begin(),
                    meshLengths.end(),
                    meshLength);
                std::copy(
                    bboxBase.begin(),
                    bboxBase.end(),
                    meshBase.begin());

                meshResolution.front() = rConfiguration["resolution"].as<std::size_t>();
                meshResolution.back() = rConfiguration["resolution"].as<std::size_t>();
            } // if discretizationType == "uniform-cartesian-grid"

            // The configuration defines a cartesian grid with potentially
            // different number of cells in each coordinate direction.
            else if (discretizationType == "cartesian-grid-2d") {
                const auto rangeConfiguration = rConfiguration["range"];
                if (rangeConfiguration.is<std::nullptr_t>()) {
                    std::copy(
                        bboxLengths.begin(),
                        bboxLengths.end(),
                        meshLengths.begin());
                    std::copy(
                        bboxBase.begin(),
                        bboxBase.end(),
                        meshBase.begin());
                } /*adaptive range*/ else {
                    for (unsigned iDimension=0u; iDimension<Dimension; ++iDimension) {
                        meshBase[iDimension] = rangeConfiguration[0][iDimension].as<double>();
                        meshLengths[iDimension] = rangeConfiguration[1][iDimension].as<double>() - meshBase[iDimension];
                    }
                }

                for (unsigned iDimension=0u; iDimension<Dimension; ++iDimension)
                    meshResolution[iDimension] = rConfiguration["resolution"][iDimension].as<std::size_t>();
            } // if discretizationType == "cartesian-grid-2d"

            else CIE_THROW(Exception, std::format(
                "unsupported geometric discretization type \"{}\"",
                discretizationType))
        CIE_END_EXCEPTION_TRACING
}


void generateMesh(
    Ref<Mesh> rMesh,
    std::span<const Scalar,Dimension> bboxBase,
    std::span<const Scalar,Dimension> bboxLengths,
    RightRef<std::vector<std::pair<
        MeshData::DomainData,
        std::vector<Scalar>>>
    > rDomainTriangles,
    std::span<const std::pair<MeshData::DomainData,Scalar>> domainMap,
    Ref<const cie::io::JSONObject> rConfiguration) {
        auto logBlock = utils::LoggerSingleton::get().newBlock("generate mesh");

        std::array<std::size_t,Dimension> meshResolution {0ul, 0ul};
        std::array<Scalar,Dimension> meshBase, meshLengths;
        parseGeometricDiscretization(
            rConfiguration["geometric"],
            bboxBase,
            bboxLengths,
            meshBase,
            meshLengths,
            meshResolution);
        std::cout << std::format(
            "mesh covers\n\tx in [{}, {}]\n\ty in [{}, {}]\n",
            meshBase.front(), meshBase.front() + meshLengths.front(),
            meshBase.back(), meshBase.back() + meshLengths.back());

        // Define an ansatz space and its derivatives.
        // In this example, every cell will use the same ansatz space.
        Ansatz ansatzSpace;

        std::array<Basis,polynomialOrder+1> basisFunctions;
        std::array<Basis::Derivative,polynomialOrder+1> basisDerivatives;

        {
            auto logBlock = utils::LoggerSingleton::get().newBlock("generate basis functions");

            enum class BasisType {
                IntegratedLegendre,
                Lagrange};
            BasisType basisType = BasisType::IntegratedLegendre;

            CIE_CHECK(
                rConfiguration["functional"]["order"].as<std::size_t>() == polynomialOrder,
                std::format(
                    "requesting a basis with a polynomial order of {}, but it is statically set to {}",
                    rConfiguration["functional"]["order"].as<std::size_t>(),
                    polynomialOrder))

            {
                const std::string basisName = rConfiguration["functional"]["type"].as<std::string>();
                if (basisName == "legendre") basisType = BasisType::IntegratedLegendre;
                else if (basisName == "lagrange") basisType = BasisType::Lagrange;
                else CIE_THROW(Exception, "unhandled basis type '" << basisName << "'")
            }

            // Construct basis functions.
            for (unsigned iBasis=0u; iBasis<polynomialOrder+1; ++iBasis) {
                // Generate a 1D polynomial serving as one of the basis functions.
                std::array<Scalar,Basis::coefficientCount> polynomialCoefficients;

                if (basisType == BasisType::IntegratedLegendre) {
                    maths::IntegratedLegendrePolynomial<Scalar> basis(iBasis);
                    CIE_CHECK(
                    basis.coefficients().size() <= polynomialCoefficients.size(),
                        "basis function " << iBasis << " is expected to have at most "
                        << polynomialCoefficients.size() << " coefficients, but has "
                        << basis.coefficients().size())
                    std::copy_n(
                        basis.coefficients().data(),
                        basis.coefficients().size(),
                        polynomialCoefficients.data());
                    std::fill_n(
                        polynomialCoefficients.data() + basis.coefficients().size(),
                        basis.coefficients().size() < polynomialCoefficients.size()
                            ? polynomialCoefficients.size() - basis.coefficients().size()
                            : 0u,
                        0.0);
                } else if (basisType == BasisType::Lagrange) {
                    std::array<Scalar,polynomialOrder + 1> nodes;
                    const Scalar nodeAngle = std::numbers::pi / polynomialOrder;
                    for (std::size_t iNode=0ul; iNode<nodes.size(); ++iNode) nodes[iNode] = std::cos(iNode * nodeAngle);
                    maths::LagrangePolynomial<Scalar> basis(
                        nodes,
                        iBasis);
                    CIE_CHECK(
                        basis.coefficients().size() <= polynomialCoefficients.size(),
                        "basis function " << iBasis << " is expected to have at most "
                        << polynomialCoefficients.size() << " coefficients, but has "
                        << basis.coefficients().size())
                    std::copy_n(
                        basis.coefficients().data(),
                        basis.coefficients().size(),
                        polynomialCoefficients.data());
                    std::fill_n(
                        polynomialCoefficients.data() + basis.coefficients().size(),
                        basis.coefficients().size() < polynomialCoefficients.size()
                            ? polynomialCoefficients.size() - basis.coefficients().size()
                            : 0u,
                        0.0);
                }

                std::cout << "basis " << iBasis << " [";
                for (auto c : polynomialCoefficients) std::cout << c << ",";
                std::cout << "],\n";
                basisFunctions[iBasis] = Basis(polynomialCoefficients);
            } // for iBasis in range(polynomialOrder + 1)

            // Construct the derivatives of all basis functions.
            for (unsigned iBasis=0u; iBasis<polynomialOrder+1; ++iBasis) {
                Ref<Basis> rBasis = basisFunctions[iBasis];
                basisDerivatives[iBasis] = rBasis.makeDerivative();
            } // for iBasis in range(polynomialOrder + 1)

            ansatzSpace = Ansatz(basisFunctions);
        }

        rMesh.data() = MeshData(
            std::move(ansatzSpace),
            std::move(rDomainTriangles),
            domainMap,
            rConfiguration);

        // Insert cells into the adjacency graph
        std::vector<Cell> cells;
        cells.reserve(std::accumulate(
            meshResolution.begin(),
            meshResolution.end(),
            1ul,
            std::multiplies<std::size_t>()));

        Size iBoundary = 0ul;
        std::array<Scalar,Dimension> edgeLengths;
        std::transform(
            meshLengths.begin(),
            meshLengths.end(),
            meshResolution.begin(),
            edgeLengths.begin(),
            std::divides<Scalar>());

        const bool discardDefaultCells = rConfiguration["geometric"]["discard-default-cells"].as<bool>();

        for (Size iCellRow : std::ranges::views::iota(0ul, meshResolution[1])) {
            for (Size iCellColumn : std::ranges::views::iota(0ul, meshResolution[0])) {
                StaticArray<StaticArray<Scalar,Dimension>,2> transformed;
                OrientedAxes<Dimension> axes;

                // Define the cell's orientation in topological and physical space.
                if (iCellColumn % 2) {
                    axes[0] = "-x";
                    transformed[0][0] = meshBase[0] + (iCellColumn + 1.0) * edgeLengths[0];
                    transformed[1][0] = meshBase[0] + iCellColumn * edgeLengths[0];
                } else {
                    axes[0] = "+x";
                    transformed[0][0] = meshBase[0] + iCellColumn * edgeLengths[0];
                    transformed[1][0] = meshBase[0] + (iCellColumn + 1.0) * edgeLengths[0];
                }

                if (iCellRow % 2) {
                    axes[1] = "-y";
                    transformed[0][1] = meshBase[1] + (iCellRow + 1.0) * edgeLengths[1];
                    transformed[1][1] = meshBase[1] + iCellRow * edgeLengths[1];
                } else {
                    axes[1] = "+y";
                    transformed[0][1] = meshBase[1] + iCellRow * edgeLengths[1];
                    transformed[1][1] = meshBase[1] + (iCellRow + 1.0) * edgeLengths[1];
                }

                SpatialTransform transform(transformed.begin(), transformed.end());
                bool isInDefaultDomain = false;
                {
                    // Ignore the cell if it lies completely inside the default domain.
                    std::array<Scalar,Dimension*intPow(2,Dimension)> corners;
                    std::vector<Scalar> buffer;
                    buffer.resize(transform.bufferSize());
                    std::size_t iCorner = 0ul;
                    ParametricSpace<Dimension,Scalar,ParametricSpaceType::Cartesian>::iterateCorners(
                        [&] (std::span<const std::uint8_t,Dimension> state) -> bool {
                            CIE_CHECK(iCorner < intPow(2,Dimension), "")
                            std::array<Scalar,Dimension> parametricCorner;
                            std::transform(
                                state.begin(),
                                state.end(),
                                parametricCorner.begin(),
                                [] (std::uint8_t c) -> Scalar {return c ? 1 : -1;});
                            transform.evaluate(
                                parametricCorner,
                                {corners.data() + iCorner * Dimension, Dimension},
                                buffer);
                            ++iCorner;
                            return true;
                        });

                    std::array<MeshData::DomainData,intPow(2,Dimension)> subdomains;
                    rMesh.data().subdomain(
                        corners,
                        subdomains);
                    isInDefaultDomain = std::none_of(
                        subdomains.begin(),
                        subdomains.end(),
                        [] (auto s) -> bool {return s;});
                }

                if (discardDefaultCells && isInDefaultDomain) continue;
                //(void)isInDefaultDomain;

                // Insert the cell into the adjacency graph (mesh) as a vertex
                const std::size_t iCell = iCellRow * meshResolution[0] + iCellColumn;
                rMesh.insert(Mesh::Vertex(
                    VertexID(iCell),
                    {} /*edges of the adjacency graph are added automatically during edge insertion*/));
                cells.emplace_back(
                    VertexID(iCell),
                    0ul,
                    1.0,
                    axes,
                    std::move(transform));

                // Insert the current cell's connections to other cells already in the mesh.
                // The rule here is that cells with lower manhattan distance from the origin
                // are sources, while those with a higher norm are targets.
                if (iCellRow) {
                    const std::size_t iSourceCell = iCell - meshResolution[0];
                    const std::size_t iTargetCell = iCell;

                    if (rMesh.findVertex(iSourceCell).has_value()) {
                        BoundaryID sharedBoundary = iCellRow % 2 ? BoundaryID("+y") : BoundaryID("-y");
                        rMesh.insert(Mesh::Edge(
                            EdgeID(iBoundary++),
                            {iSourceCell, iTargetCell},
                            {sharedBoundary}));
                    }
                } // if iCellRow

                if (iCellColumn) {
                    const std::size_t iSourceCell = iCell - 1ul;
                    const std::size_t iTargetCell = iCell;

                    if (rMesh.findVertex(iSourceCell).has_value()) {
                        BoundaryID sharedBoundary = iCellColumn % 2 ? BoundaryID("+x") : BoundaryID("-x");
                        rMesh.insert(Mesh::Edge(
                            EdgeID(iBoundary++),
                            {iSourceCell, iTargetCell},
                            {sharedBoundary}
                        ));
                    }
                } // if iCellColumn
            } // for iCellColumn in range(nodesPerDirection -1)
        } // for iCellRow in range(nodesPerDirection - 1)

        rMesh.data().setCells(std::move(cells));
        std::cout << "generated " << rMesh.vertices().size() << " cells\n";
}


} // namespace cie::fem
