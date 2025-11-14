// --- External Includes ---
#include "Eigen/Sparse"
#include "Eigen/IterativeLinearSolvers"

// --- Utility Includes ---
#include "packages/testing/inc/essentials.hpp"
#include "packages/io/inc/MatrixMarket.hpp"
#include "packages/maths/inc/Comparison.hpp"

// --- FEM Includes ---
#include "packages/graph/inc/OrientedBoundary.hpp"
#include "packages/graph/inc/OrientedAxes.hpp"
#include "packages/graph/inc/Graph.hpp"
#include "packages/graph/inc/BoundaryID.hpp"
#include "packages/graph/inc/connectivity.hpp"
#include "packages/graph/inc/Assembler.hpp"
#include "packages/maths/inc/LegendrePolynomial.hpp"
#include "packages/maths/inc/AnsatzSpace.hpp"
#include "packages/maths/inc/ScaleTranslateTransform.hpp"
#include "packages/maths/inc/AffineEmbedding.hpp"
#include "packages/maths/inc/LambdaExpression.hpp"
#include "packages/numeric/inc/GaussLegendreQuadrature.hpp"
#include "packages/numeric/inc/Quadrature.hpp"
//#include "packages/io/inc/Graphviz.hpp"
#include "packages/io/inc/GraphML.hpp"
#include "packages/io/inc/GraphML_specializations.hpp"
#include "packages/integrands/inc/LinearIsotropicStiffnessIntegrand.hpp"
#include "packages/integrands/inc/DirichletPenaltyIntegrand.hpp"
#include "packages/integrands/inc/TransformedIntegrand.hpp"

// --- GEO Includes ---
#include "packages/partitioning/inc/AABBoxNode.hpp"
#include "packages/trees/inc//ContiguousSpaceTree.hpp"

// --- STL Includes ---
#include <packages/partitioning/inc/BoxBoundable.hpp>
#include <packages/primitives/inc/Object.hpp>
#include <ranges> // ranges::iota
#include <numbers> // std::numbers::pi


namespace cie::fem {


/// @brief Number of nodes in each direction of the discretized domain.
constexpr unsigned nodesPerDirection        = 30;

constexpr unsigned polynomialOrder          = 7;

/// @brief Quadrature order used for integrands over the domain.
constexpr unsigned integrationOrder         = polynomialOrder;

/// @brief Radius of the circle Dirichlet conditions are imposed on.
constexpr double boundaryRadius             = 2.5e-1;

/// @brief Number of nodes the boundary circle is discretized by.
constexpr unsigned boundaryResolution       = 20;

/// @brief Quadrature order used for integrands over the boundary.
constexpr unsigned boundaryIntegrationOrder = polynomialOrder;

/// @brief Minimum depth of the spatial tree used to find intersections between domain cells and boundary cells.
constexpr unsigned minBoundaryTreeDepth     = 3;

/// @brief Maximum depth of the spatial tree used to find intersections between domain cells and boundary cells.
constexpr unsigned maxBoundaryTreeDepth     = 20;

/// @brief Minimum norm of a boundary segment to integrate over.
constexpr double minBoundarySegmentNorm     = 1e-12;

/// @brief Penalty parameter used for the weak imposition of Dirichlet conditions.
constexpr double weakDirichletPenalty       = 1e-5 * nodesPerDirection * nodesPerDirection;

/// @brief Number of sample points in each direction of every element to postprocess the solution on.
constexpr unsigned postprocessResolution    = 1e2;

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

    VertexID id;

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

    SpatialTransform::Inverse inverseSpatialTransform;

    CellData() noexcept = default;

    CellData(VertexID id_,
             unsigned short iAnsatz_,
             Scalar diffusivity_,
             OrientedAxes<Dimension> axes_,
             RightRef<SpatialTransform> rSpatialTransform) noexcept
        : id(id_),
          iAnsatz(iAnsatz_),
          diffusivity(diffusivity_),
          axes(axes_),
          spatialTransform(std::move(rSpatialTransform))
    {
        this->inverseSpatialTransform = spatialTransform.makeInverse();
    }

    bool at(geo::BoxBoundable<Dimension,Scalar>::Point point) const
    {
        const utils::Comparison<Scalar> comparison(1e-8, 1e-10);
        StaticArray<Scalar,Dimension> local;

        this->inverseSpatialTransform.evaluate(point, local);

        return std::all_of(
            local.begin(),
            local.end(),
            [&comparison](Scalar coordinate) {
                return comparison.less(std::abs(coordinate), static_cast<Scalar>(1))
                    || comparison.equal(std::abs(coordinate), static_cast<Scalar>(1));}
        );
    }

protected:
    void computeBoundingBoxImpl(BoundingBox& rBox) noexcept override
    {
        BoundingBox::Point opposite, localCorner, globalCorner;

        std::fill(rBox.base().begin(), rBox.base().end(), std::numeric_limits<Scalar>::max());
        std::fill(opposite.begin(), opposite.end(), std::numeric_limits<Scalar>::lowest());

        StaticArray<std::uint8_t,Dimension> state {0u, 0u};
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
        } while (cie::maths::OuterProduct<Dimension>::next(2u, state.data()));

        std::transform(opposite.begin(),
                       opposite.end(),
                       rBox.base().begin(),
                       rBox.lengths().begin(),
                       std::minus<Scalar>());
    }
}; // struct CellData


static_assert(::cie::concepts::SamplableGeometry<CellData>);


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


/// @brief Data structure unique to the triangulated, immersed boundary cells.
using BoundaryCellData = maths::AffineEmbedding<Scalar,1u,Dimension>;


/// @brief Data structure unqiue to the triangulated, immersed boundary cell corners.
using BoundaryCornerData = Scalar;


/// @brief Mesh type of the immersed, triangulated boundary.
/// @details Cell data consists of
using BoundaryMesh = Graph<
    BoundaryCellData,
    BoundaryCornerData
>;


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
        rThis.instance().inverseSpatialTransform = rThis.instance().spatialTransform.makeInverse();
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
                  Size nodesPerDirection,
                  Ref<const mp::ThreadPoolBase> rThreadPool)
{
    // Define an ansatz space and its derivatives.
    // In this example, every cell will use the same ansatz space.
    {
        Ansatz::AnsatzSet ansatzSet;
        for (unsigned iBasis=0; iBasis<polynomialOrder + 1; ++iBasis) {
            Basis::Coefficients coefficients;
            maths::IntegratedLegendrePolynomial<Scalar> legendre(iBasis);
            std::ranges::copy(legendre.coefficients(), std::back_inserter(coefficients));
            ansatzSet.emplace_back(coefficients);
        }
        rMesh.data().ansatzSpaces.emplace_back(std::move(ansatzSet), rThreadPool);
    }

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
                VertexID(iCell), // <= todo: remove duplicate id
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
    constexpr int targetLeafWidth = 5;
    constexpr int maxTreeDepth = 5;
    constexpr Scalar epsilon = 1e-3;

    geo::AABBoxNode<CellData> root;

    geo::AABBoxNode<CellData>::Point rootBase, rootLengths;
    std::fill(rootBase.begin(), rootBase.end(), -epsilon);
    std::fill(rootLengths.begin(), rootLengths.end(), 1.0 + 3.0 * epsilon);
    root = geo::AABBoxNode<CellData>(rootBase, rootLengths, nullptr);

    CIE_TEST_CASE_INIT("make bounding volume hierarchy")

    for (auto& rCell : rMesh.vertices()) {
        root.insert(&rCell.data());
    }

    root.partition(targetLeafWidth, maxTreeDepth);
    root.shrink();

    CIE_TEST_CHECK(!root.isLeaf());

    return root;
}


BoundaryMesh generateBoundaryMesh(const unsigned resolution)
{
    using Point = maths::AffineEmbedding<Scalar,1u,Dimension>::OutPoint;

    if (resolution < 3) CIE_THROW(Exception, "boundary resolution must be 3 or greater")

    BoundaryMesh boundary;

    DynamicArray<Point> corners;
    for (unsigned iSegment=0u; iSegment<resolution + 1; ++iSegment) {
        const Scalar arcParameter = iSegment * 2 * std::numbers::pi / resolution;
        corners.push_back(Point {
            boundaryRadius * std::cos(arcParameter) + 0.5,
            boundaryRadius * std::sin(arcParameter) + 0.5
        });
    } // for iSegment in range(resolution)

    for (unsigned iCorner=0u; iCorner<corners.size(); ++iCorner) {
        const StaticArray<Point,2> transformed {
            corners[iCorner],
            corners[(iCorner + 1) % corners.size()]
        };

        boundary.insert(BoundaryMesh::Vertex(
            boundary.vertices().size(),
            {},
            maths::AffineEmbedding<Scalar,1u,Dimension>(transformed)
        ));
    } // for iCorner in range(1, corners.size())

    for (unsigned iCorner=0; iCorner<corners.size(); ++iCorner) {
        boundary.insert(BoundaryMesh::Edge(
            iCorner,
            {iCorner, (iCorner + 1) % boundary.vertices().size()},
            static_cast<Scalar>(iCorner)
        ));
    } // for iCorner in range(corners.size())

    return boundary;
}


Ptr<const CellData> findContainingCell(Ref<const geo::AABBoxNode<CellData>> rBVH,
                                       Ref<const geo::AABBoxNode<CellData>::Point> rPoint)
{
    const auto pBVHNode = rBVH.find(rPoint);
    if (!pBVHNode) return nullptr;

    CIE_TEST_CHECK(pBVHNode->contains(geo::boundingBox(rPoint)));

    for (const auto pCellData : pBVHNode->contained()) {
        if (pCellData->at(rPoint)) {
            return pCellData;
        }
    }

    for (const auto pCellData : pBVHNode->intersected()) {
        if (pCellData->at(rPoint)) {
            return pCellData;
        }
    }

    CIE_TEST_CHECK(false);
    return nullptr;
}


template <class TDofMap>
void addLHSContribution(std::span<const Scalar> contribution,
                        TDofMap dofMap,
                        std::span<const int> rowExtents,
                        std::span<const int> columnIndices,
                        std::span<Scalar> entries)
{
    const unsigned localSystemSize = dofMap.size();
    for (unsigned iLocalRow=0u; iLocalRow<localSystemSize; ++iLocalRow) {
        for (unsigned iLocalColumn=0u; iLocalColumn<localSystemSize; ++iLocalColumn) {
            CIE_TEST_REQUIRE(iLocalRow < dofMap.size());
            CIE_TEST_REQUIRE(iLocalColumn < dofMap.size());

            const auto iRowBegin = rowExtents[dofMap[iLocalRow]];
            const auto iRowEnd = rowExtents[dofMap[iLocalRow] + 1];
            const auto itColumnIndex = std::lower_bound(columnIndices.begin() + iRowBegin,
                                                        columnIndices.begin() + iRowEnd,
                                                        dofMap[iLocalColumn]);
            CIE_OUT_OF_RANGE_CHECK(itColumnIndex != columnIndices.begin() + iRowEnd
                                && *itColumnIndex == dofMap[iLocalColumn]);
            const auto iEntry = std::distance(columnIndices.begin(),
                                            itColumnIndex);
            entries[iEntry] += contribution[iLocalRow * localSystemSize + iLocalColumn];
        } // for iLocalColumn in range(ansatzBuffer.size)
    } // for iLocalRow in range(ansatzBuffer.size)
}


template <class TDofMap>
void addRHSContribution(std::span<const Scalar> contribution,
                        TDofMap dofMap,
                        std::span<Scalar> rhs)
{
    for (unsigned iComponent=0u; iComponent<contribution.size(); ++iComponent) {
        const auto iRow = dofMap[iComponent];
        rhs[iRow] += contribution[iComponent];
    }
}


struct DirichletBoundary : public maths::ExpressionTraits<Scalar>
{
    using maths::ExpressionTraits<Scalar>::Span;
    using maths::ExpressionTraits<Scalar>::ConstSpan;

    void evaluate([[maybe_unused]] ConstSpan position, Span state) const noexcept
    {
        CIE_TEST_CHECK(position.size() == Dimension);
        CIE_TEST_CHECK(state.size() == 1);
        state[0] = position[0] + position[1];
    }

    unsigned size() const noexcept
    {
        return 1;
    }
}; // class DirichletBoundary


[[nodiscard]] DynamicArray<StaticArray<Scalar,2*Dimension+1>>
imposeBoundaryConditions(Ref<Mesh> rMesh,
                         Ref<const Assembler> rAssembler,
                         Ref<const geo::AABBoxNode<CellData>> rBVH,
                         std::span<const int> rowExtents,
                         std::span<const int> columnIndices,
                         std::span<Scalar> entries,
                         std::span<Scalar> rhs)
{
    DynamicArray<StaticArray<Scalar,2*Dimension+1>> boundarySegments; // {p0x, p0y, p1x, p1y, level}
    CIE_TEST_CASE_INIT("weak boundary condition imposition")

    // Load the boundary mesh.
    const auto boundary = generateBoundaryMesh(boundaryResolution);

    using TreePrimitive = geo::Cube<1,Scalar>;
    using Tree          = geo::ContiguousSpaceTree<TreePrimitive,unsigned>;

    const Quadrature<Scalar,1> lineQuadrature((GaussLegendreQuadrature<Scalar>(boundaryIntegrationOrder)));
    DynamicArray<Scalar> quadratureBuffer(  std::pow(rMesh.data().ansatzSpaces.front().size(),Dimension)
                                          + rMesh.data().ansatzSpaces.front().size());
    DynamicArray<Scalar> integrandBuffer(  rMesh.data().ansatzSpaces.front().size()
                                         + Dimension
                                         + 1);

    DirichletBoundary dirichletBoundary;
    Scalar boundaryLength = 0.0;

    for (const auto& rBoundaryCell : boundary.vertices()) {
        Tree tree(Tree::Point {-1.0}, 2.0);

        // Define a functor detecting cell boundaries.
        const auto boundaryVisitor = [&rBVH, &tree, &rBoundaryCell,
                                      &lineQuadrature, &integrandBuffer, &quadratureBuffer, &rMesh,
                                      &rowExtents, &columnIndices, &entries, &rAssembler,
                                      &rhs, &dirichletBoundary, &boundaryLength,
                                      &boundarySegments] (
            Ref<const Tree::Node> rNode,
            unsigned level) -> bool {

            if (level < minBoundaryTreeDepth) return true;
            if (maxBoundaryTreeDepth < level) return false;

            Scalar base, edgeLength;
            tree.getNodeGeometry(rNode, &base, &edgeLength);

            // Define sample points in the boundary cell's local space.
            const StaticArray<Scalar,1>
                localBase {base},
                localOpposite {base + edgeLength};

            // Transform sample points to the global coordinate space.
            geo::Traits<Dimension,Scalar>::Point globalBase, globalOpposite;
            rBoundaryCell.data().evaluate(localBase, globalBase);
            rBoundaryCell.data().evaluate(localOpposite, globalOpposite);

            // Check whether the two endpoints are in different cells.
            Ptr<const CellData> pBaseCell     = findContainingCell(rBVH, globalBase),
                                pOppositeCell = findContainingCell(rBVH, globalOpposite);

            // Integrate if both endpoints lie in the same cell.
            if (pBaseCell && pOppositeCell && pBaseCell == pOppositeCell) {
                const Scalar segmentNorm =   std::pow(globalOpposite[0] - globalBase[0], static_cast<Scalar>(2))
                                           + std::pow(globalOpposite[1] - globalBase[1], static_cast<Scalar>(2));

                if (minBoundarySegmentNorm < segmentNorm) {
                    Ref<const CellData> rCell = *pBaseCell;
                    const auto& rAnsatzSpace = rMesh.data().ansatzSpaces[rCell.iAnsatz];

                    StaticArray<maths::AffineEmbedding<Scalar,1,Dimension>::OutPoint,2> globalCorners;
                    globalCorners[0][0] = globalBase[0];
                    globalCorners[0][1] = globalBase[1];
                    globalCorners[1][0] = globalOpposite[0];
                    globalCorners[1][1] = globalOpposite[1];
                    const maths::AffineEmbedding<Scalar,1,Dimension> segmentTransform(globalCorners);

                    const auto integrand = makeTransformedIntegrand(
                        makeDirichletPenaltyIntegrand(dirichletBoundary,
                                                      /*penalty=*/weakDirichletPenalty,
                                                      rAnsatzSpace,
                                                      segmentTransform,
                                                      std::span<Scalar>(integrandBuffer)),
                        segmentTransform.makeDerivative()
                    );
                    lineQuadrature.evaluate(integrand, quadratureBuffer);

                    {
                        const auto keys = rAssembler.keys();
                        CIE_TEST_REQUIRE(std::find(keys.begin(), keys.end(), rCell.id) != keys.end());
                    }
                    const auto& rGlobalDofIndices = rAssembler[rCell.id];

                    addLHSContribution({quadratureBuffer.data(),
                                           static_cast<std::size_t>(std::pow(rAnsatzSpace.size(),Dimension))},
                                       rGlobalDofIndices,
                                       rowExtents,
                                       columnIndices,
                                       entries);

                    addRHSContribution({quadratureBuffer.data() + static_cast<std::size_t>(std::pow(rAnsatzSpace.size(),Dimension)),
                                           quadratureBuffer.data() + quadratureBuffer.size()},
                                       rGlobalDofIndices,
                                       rhs);

                    // Log debug and output info.
                    boundaryLength += std::sqrt(segmentNorm);
                    decltype(boundarySegments)::value_type segment;
                    std::copy_n(globalCorners[0].data(), Dimension, segment.data());
                    std::copy_n(globalCorners[1].data(), Dimension, segment.data() + Dimension);
                    segment.back() = level;
                    boundarySegments.push_back(segment);
                }
            } // if both endpoints lie in the same cell

            return pBaseCell != pOppositeCell && (pBaseCell || pOppositeCell);
        }; // boundaryVisitor

        // Construct a binary tree that detects intersections between
        // Cell boundaries and the current boundary cell.
        CIE_TEST_CHECK_NOTHROW(tree.scan(boundaryVisitor));
    } // for rBoundaryCell in boundary.vertices()

    CIE_TEST_CHECK(boundaryLength == Approx(boundaryRadius * 2 *std::numbers::pi).margin(5e-2));
    return boundarySegments;
}


/// @brief 2D system test.
CIE_TEST_CASE("2D", "[systemTests]")
{
    CIE_TEST_CASE_INIT("2D")

    Mesh mesh;
    mp::ThreadPoolBase threads;

    // Fill the mesh with cells and boundaries.
    {
        CIE_TEST_CASE_INIT("generate mesh")
        generateMesh(mesh, nodesPerDirection, threads);
    }

    {
        CIE_TEST_CASE_INIT("write graphml")
        io::GraphML::GraphML::Output("test_2d.graphml")(mesh);
    }

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
        DynamicArray<Scalar> integrandBuffer(std::pow(mesh.data().ansatzSpaces.front().size(), Dimension));
        DynamicArray<Scalar> productBuffer(integrandBuffer.size());

        for (Ref<const Mesh::Vertex> rCell : mesh.vertices()) {
            auto& rAnsatzDerivatives  = mesh.data().ansatzDerivatives[rCell.data().iAnsatz];
            const auto jacobian = rCell.data().spatialTransform.makeDerivative();

            const auto localIntegrand = makeTransformedIntegrand(
                LinearIsotropicStiffnessIntegrand<Ansatz::Derivative>(rCell.data().diffusivity,
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
            addLHSContribution(integrandBuffer,
                               rGlobalDofIndices,
                               rowExtents,
                               columnIndices,
                               entries);
        } // for rCell in mesh.vertices
    } // integrate

    auto bvh = makeBoundingVolumeHierarchy(mesh);
    const auto boundarySegments = imposeBoundaryConditions(
        mesh,
        assembler,
        bvh,
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

    {
        CIE_TEST_CASE_INIT("write RHS vector")

        std::cout << "rhs norm: "
                  << std::inner_product(rhs.begin(), rhs.end(), rhs.begin(), 0.0)
                  << std::endl;
        std::ofstream file("rhs.mm");
        cie::io::MatrixMarket::Output io(file);
        io(rhs.data(), rhs.size());
    }

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

    // XDMF output.
    {
        // Collect output for the bounding volume hierarchy.
        DynamicArray<StaticArray<Scalar,4 * Dimension + 1>> boundingVolumes; // {p0x, p0y, p1x, p1y, p2x, p2y, p3x, p3y}
        bvh.visit([&boundingVolumes](const auto& pBoundingVolume) -> bool {
            constexpr unsigned cornerCount = intPow(2,Dimension);
            StaticArray<geo::AABBoxNode<CellData>::Point,cornerCount> corners;
            boundingVolumes.emplace_back();
            pBoundingVolume->makeCorners(std::span<geo::AABBoxNode<CellData>::Point,cornerCount>{corners.data(),cornerCount});
            for (unsigned iCorner=0u; iCorner<cornerCount; ++iCorner) {
                std::copy_n(corners[iCorner].data(),
                            Dimension,
                            boundingVolumes.back().data() + iCorner * Dimension);
            }
            boundingVolumes.back().back() = pBoundingVolume->level();
            return true;
        });

        // Collect state samples.
        struct Sample {
            StaticArray<Scalar,Dimension>   position;
            Scalar                          state;
            unsigned                        cellID;
        }; // struct Sample
        DynamicArray<Sample> samples(postprocessResolution * postprocessResolution);

        {
            CIE_TEST_CASE_INIT("scatter postprocess")
            constexpr Scalar epsilon = 1e-10;
            constexpr Scalar postprocessDelta  = (1.0 - 2 * epsilon) / (postprocessResolution - 1);
            DynamicArray<Scalar> ansatzBuffer(mesh.data().ansatzSpaces.front().size());

            mp::ParallelFor<unsigned>(threads).firstPrivate(DynamicArray<Scalar>())(
                intPow(postprocessResolution, 2),
                [&samples, &solution, &assembler, &mesh, &bvh](const unsigned iSample,
                                                               Ref<DynamicArray<Scalar>> rAnsatzBuffer) -> void {
                    const unsigned iSampleY = iSample / postprocessResolution;
                    const unsigned iSampleX = iSampleY % postprocessResolution;
                    auto& rSample = samples[iSample];
                    rSample.position = {epsilon + iSampleX * postprocessDelta,
                                        epsilon + iSampleY * postprocessDelta};

                    // Find which cell the global point lies in.
                    Ptr<const CellData> pCellData = findContainingCell(bvh, rSample.position);

                    if (pCellData) {
                        rSample.cellID = pCellData->id;

                        // Compute sample point in the cell's local space.
                        StaticArray<Scalar,2> localSamplePoint;
                        pCellData->inverseSpatialTransform.evaluate(rSample.position, localSamplePoint);

                        // Evaluate the cell's ansatz functions at the local sample point.
                        const auto& rAnsatzSpace = mesh.data().ansatzSpaces[pCellData->iAnsatz];
                        rAnsatzBuffer.resize(rAnsatzSpace.size());
                        rAnsatzSpace.evaluate(localSamplePoint, rAnsatzBuffer);

                        // Find the entries of the cell's DoFs in the global state vector.
                        const auto& rGlobalIndices = assembler[pCellData->id];

                        // Compute state as an indirect inner product of the solution vector
                        // and the ansatz function values at the local corner coordinates.
                        rSample.state = 0;
                        for (unsigned iFunction=0u; iFunction<rGlobalIndices.size(); ++iFunction) {
                            rSample.state += solution[rGlobalIndices[iFunction]] * rAnsatzBuffer[iFunction];
                        } // for iFunction in range(rGlobalIndices.size())
                    } else {
                        // Could not find the cell that contains the current sample point.
                        rSample.state = NAN;
                    }
                });
        }

        // Write XDMF.
        std::ofstream xdmf("test_2d.xdmf");
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

                <Attribute Name="level" Center="Cell" AttributeType="Scalar">
                    <DataItem Format="XML" Dimensions=")" << boundingVolumes.size() << R"( 1">
                        )"; {
                            for (const auto& rBoundingVolume : boundingVolumes) {
                                xdmf << rBoundingVolume.back() << "\n                        ";
                            }
                        } xdmf << R"(
                    </DataItem>
                </Attribute>
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
                            StaticArray<StaticArray<Scalar,Dimension>,intPow(2,Dimension)> localCorners {
                                {-1.0, -1.0},
                                { 1.0, -1.0},
                                {-1.0,  1.0},
                                { 1.0,  1.0}
                            };

                            for (const auto& rCell : mesh.vertices()) {
                                for (const auto& rLocalPoint : localCorners) {
                                    StaticArray<Scalar,Dimension> globalPoint;
                                    rCell.data().spatialTransform.evaluate(rLocalPoint, globalPoint);
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
                            DynamicArray<Scalar> ansatzBuffer(mesh.data().ansatzSpaces.front().size());

                            for (const auto& rCell : mesh.vertices()) {
                                const auto& rGlobalIndices = assembler[rCell.id()];
                                const auto& rAnsatzSpace = mesh.data().ansatzSpaces[rCell.data().iAnsatz];
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
                                xdmf << rSample.position[0] << " " << rSample.position[1] << " 0\n                        ";
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
}


} // namespace cie::fem
