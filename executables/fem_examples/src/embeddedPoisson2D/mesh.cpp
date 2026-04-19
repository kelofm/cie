// --- Internal Includes ---
#include "embeddedPoisson2D/mesh.hpp"

// --- FEM Includes ---
#include "packages/maths/inc/LegendrePolynomial.hpp"
#include "packages/maths/inc/LagrangePolynomial.hpp"

// --- Utility Includes ---
#include "packages/logging/inc/LoggerSingleton.hpp"
#include "packages/logging/inc/LogBlock.hpp"


namespace cie::fem {


void generateMesh(
    Ref<Mesh> rMesh,
    std::span<const Scalar,2> meshBase,
    std::span<const Scalar,2> meshLengths,
    RightRef<std::vector<std::pair<
        MeshData::DomainData,
        std::vector<Scalar>>>
    > rDomainTriangles,
    std::span<const std::pair<MeshData::DomainData,Scalar>> domainMap,
    Ref<const utils::ArgParse::Results> rArguments) {
        auto logBlock = utils::LoggerSingleton::get().newBlock("generate mesh");
        const unsigned nodesPerDirection = rArguments.get<std::size_t>("r");

        {
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

                {
                    const std::string basisName = rArguments.get<std::string>("basis");
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
                domainMap);
        }

        // Insert cells into the adjacency graph
        const std::array<Scalar,2> edgeLengths {
            meshLengths.front() / (nodesPerDirection - 1),
            meshLengths.back()  / (nodesPerDirection - 1)};
        Size iBoundary = 0ul;

        for (Size iCellRow : std::ranges::views::iota(0ul, nodesPerDirection - 1)) {
            for (Size iCellColumn : std::ranges::views::iota(0ul, nodesPerDirection - 1)) {
                StaticArray<StaticArray<Scalar,Dimension>,2> transformed;
                OrientedAxes<Dimension> axes;

                // Define the cell's orientation in topological and physical space.
                if (iCellRow % 2) {
                    axes[0] = "-x";
                    transformed[0][0] = meshBase[0] + (iCellRow + 1.0) * edgeLengths[0];
                    transformed[1][0] = meshBase[0] + iCellRow * edgeLengths[0];
                } else {
                    axes[0] = "+x";
                    transformed[0][0] = meshBase[0] + iCellRow * edgeLengths[0];
                    transformed[1][0] = meshBase[0] + (iCellRow + 1.0) * edgeLengths[0];
                }

                if (iCellColumn % 2) {
                    axes[1] = "-y";
                    transformed[0][1] = meshBase[1] + (iCellColumn + 1.0) * edgeLengths[1];
                    transformed[1][1] = meshBase[1] + iCellColumn * edgeLengths[1];
                } else {
                    axes[1] = "+y";
                    transformed[0][1] = meshBase[1] + iCellColumn * edgeLengths[1];
                    transformed[1][1] = meshBase[1] + (iCellColumn + 1.0) * edgeLengths[1];
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

                if (isInDefaultDomain) continue;

                // Insert the cell into the adjacency graph (mesh) as a vertex
                const Size iCell = iCellRow * (nodesPerDirection - 1u) + iCellColumn;
                Mesh::Vertex::Data data (
                    VertexID(iCell), // <= todo: remove duplicate id
                    0ul,
                    1.0,
                    axes,
                    std::move(transform));
                rMesh.insert(Mesh::Vertex(
                    VertexID(iCell),
                    {}, ///< edges of the adjacency graph are added automatically during edge insertion
                    std::move(data)));

                // Insert the current cell's connections to other cells already in the mesh.
                // The rule here is that cells with lower manhattan distance from the origin
                // are sources, while those with a higher norm are targets.
                if (iCellRow) {
                    const Size iSourceCell = iCell - (nodesPerDirection - 1);
                    const Size iTargetCell = iCell;

                    if (rMesh.findVertex(iSourceCell).has_value()) {
                        BoundaryID sharedBoundary = iCellRow % 2 ? BoundaryID("+x") : BoundaryID("-x");
                        rMesh.insert(Mesh::Edge(
                            EdgeID(iBoundary++),
                            {iSourceCell, iTargetCell},
                            {sharedBoundary}));
                    }
                } // if iCellRow

                if (iCellColumn) {
                    const Size iSourceCell = iCell - 1ul;
                    const Size iTargetCell = iCell;

                    if (rMesh.findVertex(iSourceCell).has_value()) {
                        BoundaryID sharedBoundary = iCellColumn % 2 ? BoundaryID("+y") : BoundaryID("-y");
                        rMesh.insert(Mesh::Edge(
                            EdgeID(iBoundary++),
                            {iSourceCell, iTargetCell},
                            {sharedBoundary}
                        ));
                    }
                } // if iCellColumn
            } // for iCellColumn in range(nodesPerDirection -1)
        } // for iCellRow in range(nodesPerDirection - 1)

        std::cout << "generated " << rMesh.vertices().size() << " cells\n";
}


} // namespace cie::fem
