// --- External Includes ---
#include "Eigen/Sparse"
#include "Eigen/IterativeLinearSolvers"

// --- Utility Includes ---
#include "packages/testing/inc/essentials.hpp"
#include "packages/io/inc/MatrixMarket.hpp"

// --- FEM Includes ---
#include "packages/graph/inc/OrientedBoundary.hpp"
#include "packages/graph/inc/OrientedAxes.hpp"
#include "packages/graph/inc/Graph.hpp"
#include "packages/graph/inc/BoundaryID.hpp"
#include "packages/graph/inc/connectivity.hpp"
#include "packages/graph/inc/Assembler.hpp"
#include "packages/maths/inc/Polynomial.hpp"
#include "packages/maths/inc/AnsatzSpace.hpp"
#include "packages/maths/inc/ScaleTranslateTransform.hpp"
#include "packages/maths/inc/AffineTransform.hpp"
#include "packages/numeric/inc/GaussLegendreQuadrature.hpp"
#include "packages/numeric/inc/Quadrature.hpp"
//#include "packages/io/inc/Graphviz.hpp"
#include "packages/io/inc/GraphML.hpp"
#include "packages/io/inc/GraphML_specializations.hpp"
#include "packages/maths/inc/LinearIsotropicStiffnessIntegrand.hpp"
#include "packages/maths/inc/TransformedIntegrand.hpp"

// --- GEO Includes ---
#include "packages/partitioning/inc/AABBoxNode.hpp"

// --- STL Includes ---
#include <ranges> // ranges::iota


namespace cie::fem {


// Settings.
constexpr Size nodesPerDirection            = 5e1;
constexpr Size integrationOrder             = 3;
constexpr unsigned postprocessResolution    = 3;


/// @brief Number of spatial dimensions the problem is defined on.
constexpr unsigned Dimension = 2u;

/// @brief Floating point scalar type to use.
using Scalar = double;

/// @brief Spatial transform type mapping cells' local space to global space.
/// @details The flexibility of this transform directly defines what kind
///          of geometries the cells can represent in a boundary conforming manner.
using SpatialTransform = maths::ScaleTranslateTransform<Scalar,Dimension>;

/// @brief Cells' basis functions' type.
using Basis = maths::Polynomial<Scalar>;

/// @brief Cells' ansatz space type.
/// @details Spans the full outer product space of the basis functions.
using Ansatz = maths::AnsatzSpace<Basis,Dimension>;


/// @brief Data structure common to the entire @ref Graph "mesh".
struct MeshData
{
    /// @brief Collection of all ansatz spaces the contained cells can refer to.
    DynamicArray<Ansatz> ansatzSpaces;

    /// @brief Collection of all ansatz spaces' derivatives the contained cells can refer to.
    DynamicArray<Ansatz::Derivative> ansatzDerivatives;
}; // struct MeshData


/// @brief Data structure unique to each @ref Graph::Vertex "cell".
struct CellData : public geo::BoxBoundable<Dimension,Scalar>
{
    using geo::BoxBoundable<Dimension,Scalar>::BoundingBox;

    /// @brief Index of the cell's ansatz space in @ref MeshData::ansatzSpaces.
    unsigned short iAnsatz;

    /// @brief Local diffusivity coefficient.
    /// @details Assumed to be constant throughout the cell for simplicity.
    Scalar diffusivity;

    /// @brief Local axis orientation of the cell to enforce continuity.
    /// @details Adjacent cells must produce identical state fields
    ///          on both sides of the shared boundary. A simple way of
    ///          ensuring this is to orient the local axes of each cell
    ///          such that matching axes point in the same direction on
    ///          both sides, except the shared boundary's normal axes,
    ///          that are pointing toward each other.
    ///          @code
    ///          + ----------- +   + ----------- +
    ///          |      y      |   |      y      |
    ///          |      ^      |   |      ^      |
    ///          |      |      |   |      |      |
    ///          |      + -> x |   | x <- +      |
    ///          |             |   |             |
    ///          |             |   |             |
    ///          + ----------- +   + ----------- +
    ///          @endcode
    ///          The oriented axes define an intermediate space between
    ///          local and global space (referred to as topological space).
    ///          DoFs that don't vanish on the shared boundaries can be
    ///          merged on both sides.
    OrientedAxes<Dimension> axes;

    /// @brief Spatial transform from local to global space.
    SpatialTransform spatialTransform;

    CellData() noexcept = default;

    CellData(unsigned short iAnsatz_,
             Scalar diffusivity_,
             OrientedAxes<Dimension> axes_,
             RightRef<SpatialTransform> rSpatialTransform) noexcept
        : iAnsatz(iAnsatz_),
          diffusivity(diffusivity_),
          axes(axes_),
          spatialTransform(std::move(rSpatialTransform))
    {}

protected:
    void computeBoundingBoxImpl(BoundingBox& rBox) noexcept override
    {
        BoundingBox::Point opposite, localCorner, globalCorner;

        std::fill(rBox.base().begin(), rBox.base().end(), std::numeric_limits<Scalar>::max());
        std::fill(opposite.begin(), opposite.end(), std::numeric_limits<Scalar>::lowest());

        StaticArray<std::uint8_t,2> state {0u, 0u};
        StaticArray<Scalar,2> ordinates {-1.0, 1.0};
        do {
            // Compute the corner in local space.
            std::transform(state.begin(),
                           state.end(),
                           localCorner.begin(),
                           [&ordinates](std::uint8_t iOrdinate) {
                                return ordinates[iOrdinate];
                           });

            // Transform the corner to global space.
            this->spatialTransform.evaluate(localCorner, globalCorner);

            // Extend box definition.
            for (unsigned iDimension=0u; iDimension<Dimension; ++iDimension) {
                rBox.base()[iDimension] = std::min(rBox.base()[iDimension], globalCorner[iDimension]);
                opposite[iDimension] = std::max(opposite[iDimension], globalCorner[iDimension]);
            } // for iDimension in range(Dimension)
        } while (maths::OuterProduct<Dimension>::next(2u, state.data()));

        std::transform(opposite.begin(),
                       opposite.end(),
                       rBox.base().begin(),
                       rBox.lengths().begin(),
                       std::minus<Scalar>());
    }
}; // struct CellData


/// @brief Data structure unique to @ref Graph::Edge "boundaries" between cells.
struct BoundaryData
{
    /// @brief Boundary identifier of the shared boundary between the adjacent cells.
    BoundaryID boundary;
}; // struct BoundaryData


/// @brief Mesh type.
/// @details The mesh is modeled by an adjacency graph whose vertices are cells,
///          and edges represent the boundaries between adjacent cells. The graph
///          is not directed, but the @ref BoundaryID of the graph edge is always
///          defined on the edge's source cell.
using Mesh = Graph<CellData,BoundaryData,MeshData>;


//using BoundaryMesh = Graph<
//
//>;


/// @brief Serializer for @ref MeshData in @p GraphML format.
template <>
struct io::GraphML::Serializer<MeshData>
{
    void header(Ref<io::GraphML::XMLElement> rElement) const
    {
        // Add default value to the header.
        io::GraphML::XMLElement defaultElement = rElement.addChild("default");
        MeshData instance;
        this->operator()(defaultElement, instance);

        // Add description to the header.
        io::GraphML::XMLElement descriptionElement = rElement.addChild("desc");
        std::stringstream description;
        description << "Data structure shared by all cells and boundaries of the mesh. "
                    << "In this case, this means the ansatz spaces of the cells as "
                    << "well as their derivatives. Each cell stores an index referring "
                    << "to their own ansatz spaces.",
        descriptionElement.setValue(description.view());
    }

    void operator()(Ref<io::GraphML::XMLElement> rElement,
                    Ref<const MeshData> rInstance) const
    {
        io::GraphML::Serializer<DynamicArray<Ansatz>> subSerializer;
        subSerializer(rElement, rInstance.ansatzSpaces);
    }
}; // struct GraphML::Serializer<MeshData>

/// @brief Deserializer for @ref MeshData in @p GraphML format.
template <>
struct io::GraphML::Deserializer<MeshData>
    : public io::GraphML::DeserializerBase<MeshData>
{
    using io::GraphML::DeserializerBase<MeshData>::DeserializerBase;

    /// @brief This function is called when an element opening tag is parsed in the XML document.
    static void onElementBegin(void* pThis,
                               std::string_view elementName,
                               std::span<io::GraphML::AttributePair>)
    {
        Ref<Deserializer> rThis = *static_cast<Ptr<Deserializer>>(pThis);

        // Defer parsing to the array of ansatz spaces.
        using SubDeserializer = io::GraphML::Deserializer<DynamicArray<Ansatz>>;
        rThis.sax().push({
            SubDeserializer::make(rThis._buffer, rThis.sax(), elementName),
            SubDeserializer::onElementBegin,
            SubDeserializer::onText,
            SubDeserializer::onElementEnd
        });
    }

    /// @brief This function is called when text block is parsed in the XML document.
    static void onText(void*,
                       std::string_view)
    {
        // No text data is expected for this class.
        CIE_THROW(Exception, "Unexpected text block while parsing mesh data.")
    }

    /// @brief This function is called when an element closing tag is parsed in the XML document.
    static void onElementEnd(void* pThis,
                             std::string_view elementName)
    {
        Ref<Deserializer> rThis = *static_cast<Ptr<Deserializer>>(pThis);

        // Move the parsed ansatz spaces from the buffer to the mesh data instance.
        rThis.instance().ansatzSpaces = std::move(rThis._buffer);

        // Build derivatives from the parsed ansatz spaces.
        std::transform(rThis.instance().ansatzSpaces.begin(),
                       rThis.instance().ansatzSpaces.end(),
                       std::back_inserter(rThis.instance().ansatzDerivatives),
                       [] (Ref<const Ansatz> rAnsatz) -> Ansatz::Derivative {
                            return rAnsatz.makeDerivative();
                       });

        // The parser's job is done => destroy it.
        rThis.template release<Deserializer>(&rThis, elementName);
    }

private:
    DynamicArray<Ansatz> _buffer;
}; // struct GraphML::Deserializer<MeshData>


template <>
struct io::GraphML::Serializer<CellData>
{
    void header(Ref<io::GraphML::XMLElement> rElement) const
    {
        io::GraphML::XMLElement defaultElement = rElement.addChild("default");
        CellData instance;
        this->operator()(defaultElement, instance);
    }

    void operator()(Ref<io::GraphML::XMLElement> rElement,
                    Ref<const CellData> rInstance) const
    {
        rElement.addAttribute("iAnsatz", std::to_string(rInstance.iAnsatz));
        rElement.addAttribute("diffusivity", std::to_string(rInstance.diffusivity));

        std::stringstream buffer;
        buffer << rInstance.axes;
        rElement.addAttribute("axes", buffer.view());

        io::GraphML::Serializer<SpatialTransform>()(rElement, rInstance.spatialTransform);
    }
}; // struct GraphML::Serializer<CellData>


template <>
struct io::GraphML::Deserializer<CellData>
    : public io::GraphML::DeserializerBase<CellData>
{
    using io::GraphML::DeserializerBase<CellData>::DeserializerBase;

    /// @brief This function is called when an element opening tag is parsed in the XML document.
    static void onElementBegin(void* pThis,
                               std::string_view elementName,
                               std::span<io::GraphML::AttributePair> attributes)
    {
        Ref<Deserializer> rThis = *static_cast<Ptr<Deserializer>>(pThis);

        {
            const auto it = std::find_if(attributes.begin(),
                                         attributes.end(),
                                         [] (const auto pair) {
                                         return pair.first == "iAnsatz";
                                         });
            char* pEnd;
            rThis.instance().iAnsatz = std::strtoul(it->second.data(), &pEnd, 10);
        }

        {
            const auto it = std::find_if(attributes.begin(),
                                         attributes.end(),
                                         [] (const auto pair) {
                                         return pair.first == "diffusivity";
                                         });
            char* pEnd;
            rThis.instance().diffusivity = std::strtod(it->second.data(), &pEnd);
        }

        {
            const auto it = std::find_if(attributes.begin(),
                                         attributes.end(),
                                         [] (const auto pair) {
                                         return pair.first == "axes";
                                         });
            rThis.instance().axes = OrientedAxes<Dimension>(it->second);
        }

        {
            using SubSerializer = io::GraphML::Deserializer<SpatialTransform>;
            rThis.sax().push({
                SubSerializer::make(rThis.instance().spatialTransform, rThis.sax(), elementName),
                SubSerializer::onElementBegin,
                SubSerializer::onText,
                SubSerializer::onElementEnd
            });
        }

        CIE_THROW(NotImplementedException, "")
    }

    static void onText(Ptr<void>, std::string_view) noexcept {}

    static void onElementEnd(Ptr<void> pThis,
                             std::string_view elementName)
    {
        Ref<Deserializer> rThis = *static_cast<Ptr<Deserializer>>(pThis);
        rThis.template release<Deserializer>(&rThis, elementName);
    }

}; // struct GraphML::Deserializer<CellData>


template <>
struct io::GraphML::Serializer<BoundaryData>
{
    void header(Ref<io::GraphML::XMLElement> rElement) const
    {
        io::GraphML::XMLElement defaultElement = rElement.addChild("default");
        BoundaryData instance;
        this->operator()(defaultElement, instance);
    }

    void operator()(Ref<io::GraphML::XMLElement> rElement,
                    Ref<const BoundaryData> rInstance) const
    {
        std::ostringstream buffer;
        buffer << rInstance.boundary;
        rElement.addAttribute("bnd", buffer.view());
    }
}; // struct GraphML::Serializer<BoundaryData>


template <>
struct io::GraphML::Deserializer<BoundaryData>
    : public io::GraphML::DeserializerBase<BoundaryData>
{
    using io::GraphML::DeserializerBase<BoundaryData>::DeserializerBase;

    /// @brief This function is called when an element opening tag is parsed in the XML document.
    static void onElementBegin(void* pThis,
                               std::string_view,
                               std::span<io::GraphML::AttributePair> attributes)
    {
        Ref<Deserializer> rThis = *static_cast<Ptr<Deserializer>>(pThis);

        {
            const auto it = std::find_if(attributes.begin(),
                                         attributes.end(),
                                         [] (const auto pair) {
                                            return pair.first == "bnd";
                                         });
            rThis.instance().boundary = BoundaryID(it->second.data());
        }
    }

    static void onText(Ptr<void>, std::string_view) noexcept
    {
        CIE_THROW(
            Exception,
            "Unexpected text block while parsing boundary data."
        )
    }

    static void onElementEnd(Ptr<void> pThis,
                             std::string_view elementName)
    {
        Ref<Deserializer> rThis = *static_cast<Ptr<Deserializer>>(pThis);
        rThis.template release<Deserializer>(&rThis, elementName);
    }

}; // struct GraphML::Deserializer<BoundaryData>


/// @brief Generate cells and boundaries for the example problem.
/// @details Mesh:
///          @code
///            [u(0,1)=2]                                     [u(1,1)=3]
///                     +---------+                 +---------+
///                     | (m-1)n+1|                 |    mn   |
///                     +---------+                 +---------+
///                       .
///                       .
///                       .
///                     +---------+
///                     |   n+1   |
///                     +---------+
///                     +---------+---------+       +---------+
///                     |    1    |    2    |  ...  |    n    |
///                     +---------+---------+       +---------+
///            [u(0,0)=0]                                     [u(1,0)=1]
///          @endcode
void generateMesh(Ref<Mesh> rMesh,
                  Size nodesPerDirection)
{
    // Define an ansatz space and its derivatives.
    // In this example, every cell will use the same ansatz space.
    rMesh.data().ansatzSpaces.emplace_back(Ansatz::AnsatzSet {
         Basis({ 0.5,  0.5      })
        ,Basis({ 0.5, -0.5      })
        ,Basis({ 1.0,  0.0, -1.0})
    });

    rMesh.data().ansatzDerivatives.emplace_back(
        rMesh.data().ansatzSpaces.front().makeDerivative()
    );

    // Insert cells into the adjacency graph
    const Scalar edgeLength = 1.0 / (nodesPerDirection - 1);
    Size iBoundary = 0ul;

    for (Size iCellRow : std::ranges::views::iota(0ul, nodesPerDirection - 1)) {
        for (Size iCellColumn : std::ranges::views::iota(0ul, nodesPerDirection - 1)) {
            StaticArray<StaticArray<Scalar,Dimension>,2> transformed;
            OrientedAxes<Dimension> axes;

            // Define the cell's orientation in topological and physical space.
            if (iCellRow % 2) {
                axes[0] = "-x";
                transformed[0][0] = (iCellRow + 1.0) * edgeLength;
                transformed[1][0] = iCellRow * edgeLength;
            } else {
                axes[0] = "+x";
                transformed[0][0] = iCellRow * edgeLength;
                transformed[1][0] = (iCellRow + 1.0) * edgeLength;
            }

            if (iCellColumn % 2) {
                axes[1] = "-y";
                transformed[0][1] = (iCellColumn + 1.0) * edgeLength;
                transformed[1][1] = iCellColumn * edgeLength;
            } else {
                axes[1] = "+y";
                transformed[0][1] = iCellColumn * edgeLength;
                transformed[1][1] = (iCellColumn + 1.0) * edgeLength;
            }

            // Insert the cell into the adjacency graph (mesh) as a vertex
            const Size iCell = iCellRow * (nodesPerDirection - 1u) + iCellColumn;
            Mesh::Vertex::Data data (
                    0u,   // <= All cells share the same ansatz space in this example.
                    1.0,
                    axes,
                    SpatialTransform(transformed.begin(), transformed.end())
                );
            rMesh.insert(Mesh::Vertex(
                VertexID(iCell),
                {}, ///< edges of the adjacency graph are added automatically during edge insertion
                std::move(data)
            ));

            // Insert the current cell's connections to other cells already in the mesh.
            // The rule here is that cells with lower manhattan distance from the origin
            // are sources, while those with a higher norm are targets.
            if (iCellRow) {
                const Size iSourceCell = iCell - (nodesPerDirection - 1);
                const Size iTargetCell = iCell;
                BoundaryID sharedBoundary = iCellRow % 2 ? BoundaryID("+x") : BoundaryID("-x");
                rMesh.insert(Mesh::Edge(
                    EdgeID(iBoundary++),
                    {iSourceCell, iTargetCell},
                    {sharedBoundary}
                ));
            } // if iCellRow

            if (iCellColumn) {
                const Size iSourceCell = iCell - 1ul;
                const Size iTargetCell = iCell;
                BoundaryID sharedBoundary = iCellColumn % 2 ? BoundaryID("+y") : BoundaryID("-y");
                rMesh.insert(Mesh::Edge(
                    EdgeID(iBoundary++),
                    {iSourceCell, iTargetCell},
                    {sharedBoundary}
                ));
            } // if iCellColumn
        } // for iCellColumn in range(nodesPerDirection -1)
    } // for iCellRow in range(nodesPerDirection - 1)
}


geo::AABBoxNode<CellData> makeBoundingVolumeHierarchy(Ref<Mesh> rMesh)
{
    geo::AABBoxNode<CellData> root;

    geo::AABBoxNode<CellData>::Point rootBase, rootLengths;
    std::fill(rootBase.begin(), rootBase.end(), 0.0);
    std::fill(rootLengths.begin(), rootLengths.end(), 1.0);
    root = geo::AABBoxNode<CellData>(rootBase, rootLengths, nullptr);

    CIE_TEST_CASE_INIT("make bounding volume hierarchy")
    constexpr int targetLeafWidth = 10;
    constexpr int maxTreeDepth = 2;

    for (auto& rCell : rMesh.vertices()) {
        root.insert(&rCell.data());
    }

    root.partition(targetLeafWidth, maxTreeDepth);

    return root;
}


void imposeBoundaryConditions(Ref<const Mesh> rMesh,
                              Ref<const Assembler> rAssembler,
                              std::span<const int> rowExtents,
                              std::span<const int> columnIndices,
                              [[maybe_unused]] std::span<Scalar> entries,
                              std::span<Scalar> rhs)
{
    // Find DoFs to constrain:
    // 0) u(0, 0) = 0
    // 1) u(1, 0) = 1
    // 2) u(0, 1) = 2
    // 3) u(1, 1) = 3
    CIE_TEST_CASE_INIT("apply boundary conditions")

    StaticArray<std::optional<std::size_t>,4> iConstrainedDofs;
    utils::Comparison<Scalar> comparison(1e-8, 1e-6);
    const int rowCount = rowExtents.size() - 1;

    for (const auto& rCell : rMesh.vertices()) {
        // The DoFs we're looking for are nodal DoFs, meaning we need to find
        // a cell that has a vertex at a corner, and find the DoF within whose
        // basis function evaluates to 1 at that location.
        std::optional<std::size_t> maybeLocalCornerIndex;
        std::optional<std::size_t> maybeGlobalCornerIndex;

        for (const unsigned iLocalCorner : std::ranges::views::iota(0u, intPow(2u, Dimension))) {
            StaticArray<Scalar,Dimension> localCoordinates, globalCoordinates;
            localCoordinates[0] = (iLocalCorner & 1u) ? 1.0 : -1.0;
            localCoordinates[1] = (iLocalCorner & 2u) ? 1.0 : -1.0;
            rCell.data().spatialTransform.evaluate(localCoordinates, globalCoordinates);

            // Loop over global corners and see whether any coincide with the transformed point.
            for (const unsigned iGlobalCorner : std::ranges::views::iota(0u, intPow(2u, Dimension))) {
                StaticArray<Scalar,Dimension> referenceCoordinates;
                referenceCoordinates[0] = Scalar(iGlobalCorner & 1u);
                referenceCoordinates[1] = Scalar((iGlobalCorner & 2u) >> 1);
                if (comparison.equal(globalCoordinates[0], referenceCoordinates[0]) && comparison.equal(globalCoordinates[1], referenceCoordinates[1])) {
                    maybeLocalCornerIndex = iLocalCorner;
                    maybeGlobalCornerIndex = iGlobalCorner;
                    break;
                }
            } //for iGlobalCorner in range(4)

            if (maybeLocalCornerIndex.has_value()) {
                break;
            }
        } // for iLocalCorner in range(4)

        // Found a cell that has a vertex at a corner.
        if (maybeLocalCornerIndex.has_value()) {
            const auto& rAnsatzSpace = rMesh.data().ansatzSpaces[rCell.data().iAnsatz];

            StaticArray<Scalar,Dimension> localCoordinates;
            localCoordinates[0] = (maybeLocalCornerIndex.value() & 1u) ? 1.0 : -1.0;
            localCoordinates[1] = (maybeLocalCornerIndex.value() & 2u) ? 1.0 : -1.0;
            DynamicArray<Scalar> ansatzValues(rAnsatzSpace.size());
            rAnsatzSpace.evaluate(localCoordinates, ansatzValues);
            for (unsigned iAnsatz=0u; iAnsatz<ansatzValues.size(); ++iAnsatz) {
                if (comparison.equal(ansatzValues[iAnsatz], 1.0)) {
                    iConstrainedDofs[maybeGlobalCornerIndex.value()] = rAssembler[rCell.id()][iAnsatz];
                    break;
                }
            }
        }
    } // for rCell in rMesh.vertices

    // Barbaric imposition of dirichlet conditions.
    for (const auto maybeDofIndex : iConstrainedDofs) {
        CIE_TEST_REQUIRE(maybeDofIndex.has_value());
    }

    DynamicArray<std::pair<std::size_t,Scalar>> dirichletConditions;
    for (unsigned iDof=0u; iDof<iConstrainedDofs.size(); ++iDof) {
        dirichletConditions.emplace_back(iConstrainedDofs[iDof].value(), Scalar(iDof));
    }

    for (const auto& [iDof, value] : dirichletConditions) {
        for (std::size_t iRow=0ul; iRow<static_cast<std::size_t>(rowCount); ++iRow) {
            const std::size_t iEntryBegin = rowExtents[iRow];
            const std::size_t iEntryEnd   = rowExtents[iRow + 1];
            for (std::size_t iEntry=iEntryBegin; iEntry<iEntryEnd; ++iEntry) {
                const std::size_t iColumn = columnIndices[iEntry];
                if (iRow == iDof) {
                    if (iRow == iColumn) entries[iEntry] = 1.0;
                    else entries[iEntry] = 0.0;
                } else if (iColumn == iDof) {
                    rhs[iRow] -= value * entries[iEntry];
                    entries[iEntry] = 0.0;
                }
            } // for iEntry in range(iEntryBegin, iEntryEnd)
        } // for iRow in range(rowCount)
        rhs[iDof] = value;
    }
}


/// @brief 2D system test.
CIE_TEST_CASE("2D", "[systemTests]")
{
    CIE_TEST_CASE_INIT("2D")

    Mesh mesh;

    // Fill the mesh with cells and boundaries.
    {
        CIE_TEST_CASE_INIT("generate mesh")
        generateMesh(mesh, nodesPerDirection);
    }

    {
        CIE_TEST_CASE_INIT("write graphml")
        io::GraphML::GraphML::Output("test_2d.graphml")(mesh);
    }

    // Create a spatial search.
    auto bvh = makeBoundingVolumeHierarchy(mesh);
    CIE_TEST_CHECK(bvh.containedObjects().size() == mesh.vertices().size());

    // Find ansatz functions that coincide on opposite boundaries.
    // In adjacent cells, these ansatz functions will have to map
    // to the same DoF in the assembled system.
    const StaticArray<Scalar,5> samples {-1.0, -0.5, 0.0, 0.5, 1.0};
    const auto ansatzMap = makeAnsatzMap(mesh.data().ansatzSpaces.front(),
                                         samples,
                                         utils::Comparison<Scalar>(/*absoluteTolerance =*/ 1e-8,
                                                                   /*relativeTolerance =*/ 1e-6));

    // Build a factory detects the topology of the mesh and issues indices to DoFs.
    Assembler assembler;

    {
        CIE_TEST_CASE_INIT("parse mesh topology")
        assembler.addGraph(mesh,
                           [&mesh]([[maybe_unused]] Ref<const Mesh::Vertex> rVertex) -> std::size_t {
                                const Ansatz& rAnsatz = mesh.data().ansatzSpaces[rVertex.data().iAnsatz];
                                return rAnsatz.size();
                           },
                           [&ansatzMap, &mesh](Ref<const Mesh::Edge> rEdge, Assembler::DoFPairIterator it) {
                                const auto sourceAxes = mesh.find(rEdge.source()).value().data().axes;
                                const auto targetAxes = mesh.find(rEdge.target()).value().data().axes;
                                ansatzMap.getPairs(OrientedBoundary<Dimension>(sourceAxes, rEdge.data().boundary),
                                                   OrientedBoundary<Dimension>(targetAxes, rEdge.data().boundary),
                                                   it);
                           });
    } // parse mesh topology

    // Create empty CSR matrix
    int rowCount, columnCount;
    DynamicArray<int> rowExtents, columnIndices;
    DynamicArray<double> entries;
    {
        CIE_TEST_CASE_INIT("compute sparsity pattern")
        assembler.makeCSRMatrix(rowCount, columnCount, rowExtents, columnIndices, entries);
    }
    DynamicArray<Scalar> rhs(rowCount, 0.0);

    // Compute element contributions and assemble them into the matrix
    {
        CIE_TEST_CASE_INIT("integrate")

        const Quadrature<Scalar,Dimension> quadrature((GaussLegendreQuadrature<Scalar>(integrationOrder)));
        DynamicArray<Scalar> derivativeBuffer(mesh.data().ansatzDerivatives.front().size());
        DynamicArray<Scalar> integrandBuffer(std::pow(mesh.data().ansatzSpaces.front().size(), 2ul));
        DynamicArray<Scalar> productBuffer(integrandBuffer.size());

        for (Ref<const Mesh::Vertex> rCell : mesh.vertices()) {
            const auto& rAnsatzSpace        = mesh.data().ansatzSpaces[rCell.data().iAnsatz];
            auto& rAnsatzDerivatives  = mesh.data().ansatzDerivatives[rCell.data().iAnsatz];
            const auto jacobian = rCell.data().spatialTransform.makeInverse().makeDerivative();

            const auto localIntegrand = maths::makeTransformedIntegrand(
                maths::LinearIsotropicStiffnessIntegrand<Ansatz::Derivative>(rCell.data().diffusivity,
                                                                             rAnsatzDerivatives,
                                                                             {derivativeBuffer.data(), derivativeBuffer.size()}),
                jacobian
            );
            quadrature.evaluate(localIntegrand, integrandBuffer);

            {
                const auto keys = assembler.keys();
                CIE_TEST_REQUIRE(std::find(keys.begin(), keys.end(), rCell.id()) != keys.end());
            }
            const auto& rGlobalDofIndices = assembler[rCell.id()];
            const unsigned localSystemSize = rAnsatzSpace.size();

            for (unsigned iLocalRow=0u; iLocalRow<localSystemSize; ++iLocalRow) {
                for (unsigned iLocalColumn=0u; iLocalColumn<localSystemSize; ++iLocalColumn) {
                    CIE_TEST_REQUIRE(iLocalRow < rGlobalDofIndices.size());
                    CIE_TEST_REQUIRE(iLocalColumn < rGlobalDofIndices.size());

                    const auto iRowBegin = rowExtents[rGlobalDofIndices[iLocalRow]];
                    const auto iRowEnd = rowExtents[rGlobalDofIndices[iLocalRow] + 1];
                    const auto itColumnIndex = std::lower_bound(columnIndices.begin() + iRowBegin,
                                                                columnIndices.begin() + iRowEnd,
                                                                rGlobalDofIndices[iLocalColumn]);
                    CIE_OUT_OF_RANGE_CHECK(itColumnIndex != columnIndices.begin() + iRowEnd
                                           && *itColumnIndex == rGlobalDofIndices[iLocalColumn]);
                    const auto iEntry = std::distance(columnIndices.begin(),
                                                      itColumnIndex);
                    entries[iEntry] += integrandBuffer[iLocalRow * localSystemSize + iLocalColumn];
                } // for iLocalColumn in range(ansatzBuffer.size)
            } // for iLocalRow in range(ansatzBuffer.size)
        } // for rCell in mesh.vertices
    } // integrate

    imposeBoundaryConditions(mesh,
                             assembler,
                             rowExtents,
                             columnIndices,
                             entries,
                             rhs);

    // Matrix market output.
    {
        CIE_TEST_CASE_INIT("write LHS matrix")
        std::ofstream file("lhs.mm");
        cie::io::MatrixMarket::Output io(file);
        io(rowCount,
           columnCount,
           entries.size(),
           rowExtents.data(),
           columnIndices.data(),
           entries.data());
    } // matrix market output

    // Solve the linear system.
    DynamicArray<Scalar> solution(rhs.size());
    {
        CIE_TEST_CASE_INIT("solve")
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

        solver.compute(lhsAdaptor);
        solutionAdaptor = solver.solve(rhsAdaptor);

        CIE_TEST_CHECK(solver.info() == Eigen::ComputationInfo::Success);
        std::cout << solver.iterations() << " iterations "
                  << solver.error()      << " residual\n";
    } // solve

//    {
//        std::ofstream file("mesh.gv");
//        io::Graphviz::Output io(file, io::Graphviz::Settings {.directed = false});
//        io(
//            mesh,
//            [&assembler](const Mesh::Vertex& rVertex) -> std::string {
//                std::stringstream label;
//                label << rVertex.id() << ":";
//                for (const auto& rItem : assembler[rVertex.id()]) {
//                    label << rItem << ',';
//                }
//                return label.str();
//            },
//            [](const Mesh::Edge& rEdge) -> std::string {
//                std::stringstream label;
//                label << rEdge.data().first << ":" << rEdge.data().second;
//                return label.str();
//            }
//        );
//    } // write GraphViz

    // Postprocessing
    {
        CIE_TEST_CASE_INIT("scatter postprocess")

        DynamicArray<std::pair<StaticArray<Scalar,Dimension>,Scalar>> solutionSamples;
        solutionSamples.reserve(postprocessResolution * intPow(nodesPerDirection - 1u, Dimension));

        constexpr Scalar postprocessDelta  = 2.0 / (postprocessResolution + 1.0);
        constexpr Scalar postprocessOffset = -1.0 + 2.0 / (postprocessResolution + 1.0);

        DynamicArray<Scalar> ansatzBuffer(mesh.data().ansatzSpaces.size());
        DynamicArray<StaticArray<Scalar,Dimension>> localSamplePoints;

        {
            StaticArray<unsigned,Dimension> samplePointState;
            std::fill(samplePointState.begin(), samplePointState.end(), 0u);
            do {
                StaticArray<Scalar,Dimension> localSamplePoint;
                for (unsigned iDimension=0u; iDimension<Dimension; ++iDimension) {
                    localSamplePoint[iDimension] = postprocessOffset + samplePointState[iDimension] * postprocessDelta;
                }
                localSamplePoints.emplace_back(localSamplePoint);
            } while (maths::OuterProduct<Dimension>::next(postprocessResolution, samplePointState.data()));
        }

        for (const auto& rCell : mesh.vertices()) {
            const auto& rGlobalIndices = assembler[rCell.id()];
            const auto& rAnsatzSpace = mesh.data().ansatzSpaces[rCell.data().iAnsatz];

            for (const auto& localCoordinates : localSamplePoints) {
                StaticArray<Scalar,Dimension> globalSamplePoint;
                rCell.data().spatialTransform.evaluate(localCoordinates, globalSamplePoint);
                solutionSamples.emplace_back(globalSamplePoint, 0.0);
                ansatzBuffer.resize(rAnsatzSpace.size());
                rAnsatzSpace.evaluate(localCoordinates, ansatzBuffer);

                for (unsigned iFunction=0u; iFunction<ansatzBuffer.size(); ++iFunction) {
                    solutionSamples.back().second += solution[rGlobalIndices[iFunction]] * ansatzBuffer[iFunction];
                }
            } // for localCoordinates in localSamplePoints
        } // for rCell in mesh.vertices

        const std::string fileName = "test_2d_solution.csv";
        std::cout << "write scatter samples to " << fileName << std::endl;
        std::ofstream file(fileName);
        file << "cells-per-direction,"       << "postprocess-resolution\n"
             << nodesPerDirection - 1 << "," << postprocessResolution << "\n"
             << "x,y,state\n";
        for (const auto& [rSamplePoint, rValue] : solutionSamples) {
            for (const auto coordinate : rSamplePoint) file << coordinate << ',';
            file << rValue << '\n';
        } // for samplePoint, rValue in solutionSamples
    } // scatter postprocess

    // STL surface.
    if (Dimension == 2 && 2 <= postprocessResolution) {
        CIE_TEST_CASE_INIT("STL postprocess")

        // Compute the total number of triangles.
        const std::uint32_t triangleCount = std::pow(nodesPerDirection - 1, 2ul) * std::pow(2 * (postprocessResolution - 1), 2ul);
        constexpr Scalar postprocessDelta  = 2.0 / (postprocessResolution - 1.0);

        std::ofstream file("test_2d_solution.stl", std::ios::binary);

        // Write STL header.
        {
            std::array<std::uint8_t,80> buffer;
            file.write(reinterpret_cast<const char*>(buffer.data()), 80);
            file.write(reinterpret_cast<const char*>(&triangleCount), sizeof(triangleCount));
        }

        DynamicArray<StaticArray<Scalar,Dimension+1>> solutionSamples;
        DynamicArray<Scalar> ansatzBuffer(mesh.data().ansatzSpaces.size());
        DynamicArray<StaticArray<Scalar,Dimension>> localSamples;

        {
            StaticArray<unsigned,Dimension> samplePointState;
            std::fill(samplePointState.begin(), samplePointState.end(), 0u);
            do {
                StaticArray<Scalar,Dimension> localSamplePoint;
                for (unsigned iDimension=0u; iDimension<Dimension; ++iDimension) {
                    localSamplePoint[iDimension] = -1.0 + samplePointState[iDimension] * postprocessDelta;
                }
                localSamples.emplace_back(localSamplePoint);
            } while (maths::OuterProduct<Dimension>::next(postprocessResolution, samplePointState.data()));
        }

        for (const auto& rCell : mesh.vertices()) {
            solutionSamples.clear();
            const auto& rGlobalIndices = assembler[rCell.id()];
            const auto& rAnsatzSpace = mesh.data().ansatzSpaces[rCell.data().iAnsatz];

            for (const auto& localCoordinates : localSamples) {
                solutionSamples.emplace_back();

                // Compute coordinates in global space.
                rCell.data().spatialTransform.evaluate(localCoordinates, solutionSamples.back());
                solutionSamples.back().back() = static_cast<Scalar>(0);

                // Compute ansatz function values at the sample points.
                ansatzBuffer.resize(rAnsatzSpace.size());
                rAnsatzSpace.evaluate(localCoordinates, ansatzBuffer);

                // Compute state as an indirect inner product of the solution vector
                // and the ansatz function values at the local sample point coordinates.
                for (unsigned iFunction=0u; iFunction<ansatzBuffer.size(); ++iFunction) {
                    solutionSamples.back().back() += solution[rGlobalIndices[iFunction]] * ansatzBuffer[iFunction];
                }
            } // for localCoordinates in localSamples

            for (unsigned iSampleX=1u; iSampleX<postprocessResolution; ++iSampleX) {
                for (unsigned iSampleY=1u; iSampleY<postprocessResolution; ++iSampleY) {
                    // Each triangle is represented by its normal, and 3 vertices.
                    StaticArray<float,12> stlTriangle;

                    // Collect and write the first triangle.
                    {
                        StaticArray<std::size_t,3> sampleIndices {
                            iSampleY * postprocessResolution + iSampleX,
                            iSampleY * postprocessResolution + iSampleX - 1,
                            (iSampleY - 1) * postprocessResolution + iSampleX - 1
                        };

                        for (unsigned iVertex=0u; iVertex<sampleIndices.size(); ++iVertex) {
                            const auto& rSource = solutionSamples[sampleIndices[iVertex]];
                            std::span<float> target(stlTriangle.data() + (iVertex + 1) * 3, 3);
                            std::copy(rSource.begin(), rSource.end(), target.begin());
                        }

                        StaticArray<StaticArray<Scalar,3>,2> edges {
                            StaticArray<Scalar,3> {stlTriangle[ 6] - stlTriangle[3],
                                                   stlTriangle[ 7] - stlTriangle[4],
                                                   stlTriangle[ 8] - stlTriangle[5]},
                            StaticArray<Scalar,3> {stlTriangle[ 9] - stlTriangle[3],
                                                   stlTriangle[10] - stlTriangle[4],
                                                   stlTriangle[11] - stlTriangle[5]},
                        };

                        stlTriangle[0] = edges[0][1] * edges[1][2] - edges[0][2] * edges[1][1];
                        stlTriangle[1] = edges[0][2] * edges[1][0] - edges[0][0] * edges[1][2];
                        stlTriangle[2] = edges[0][0] * edges[1][1] - edges[0][1] * edges[1][0];

                        file.write(reinterpret_cast<char*>(stlTriangle.data()),
                                   sizeof(float) * stlTriangle.size());

                        std::uint16_t attribute = 0;
                        file.write(reinterpret_cast<char*>(&attribute), sizeof(attribute));
                    }

                    // Collect and write the second triangle.
                    {
                        StaticArray<std::size_t,3> sampleIndices {
                            iSampleY * postprocessResolution + iSampleX,
                            (iSampleY - 1) * postprocessResolution + iSampleX - 1,
                            (iSampleY - 1) * postprocessResolution + iSampleX,
                        };

                        for (unsigned iVertex=0u; iVertex<sampleIndices.size(); ++iVertex) {
                            const auto& rSource = solutionSamples[sampleIndices[iVertex]];
                            std::span<float> target(stlTriangle.data() + (iVertex + 1) * 3, 3);
                            std::copy(rSource.begin(), rSource.end(), target.begin());
                        }

                        StaticArray<StaticArray<Scalar,3>,2> edges {
                            StaticArray<Scalar,3> {stlTriangle[ 6] - stlTriangle[3],
                                                   stlTriangle[ 7] - stlTriangle[4],
                                                   stlTriangle[ 8] - stlTriangle[5]},
                            StaticArray<Scalar,3> {stlTriangle[ 9] - stlTriangle[3],
                                                   stlTriangle[10] - stlTriangle[4],
                                                   stlTriangle[11] - stlTriangle[5]},
                        };

                        stlTriangle[0] = edges[0][1] * edges[1][2] - edges[0][2] * edges[1][1];
                        stlTriangle[1] = edges[0][2] * edges[1][0] - edges[0][0] * edges[1][2];
                        stlTriangle[2] = edges[0][0] * edges[1][1] - edges[0][1] * edges[1][0];

                        file.write(reinterpret_cast<char*>(stlTriangle.data()),
                                   sizeof(float) * stlTriangle.size());

                        std::uint16_t attribute = 0;
                        file.write(reinterpret_cast<char*>(&attribute), sizeof(attribute));
                    }

                    // Collect the second triangle.
                } // for iSampleY in range(1, postprocessResolution)
            } // for iSampleX in range(1, postprocessResolution)
        } // for rCell in mesh.vertices
    } // STL postprocess
}


} // namespace cie::fem
