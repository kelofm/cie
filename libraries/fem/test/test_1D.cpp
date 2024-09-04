// --- Utility Includes ---
#include "packages/testing/inc/essentials.hpp"
#include "packages/stl_extension/inc/SymmetricPair.hpp"
#include "packages/io/inc/MatrixMarket.hpp"

// --- FEM Includes ---
#include "packages/graph/inc/Graph.hpp"
#include "packages/graph/inc/BoundaryID.hpp"
#include "packages/graph/inc/connectivity.hpp"
#include "packages/graph/inc/Assembler.hpp"
#include "packages/maths/inc/Polynomial.hpp"
#include "packages/maths/inc/AnsatzSpace.hpp"
#include "packages/utilities/inc/AttributeContainer.hpp"

// --- STL Includes ---
#include <ranges> // ranges::iota


namespace cie::fem {


using Scalar = double;
using CellID = Size;
using NodeID = Size;
using DoFID  = Size;


template <class TAnsatzSpace>
struct CellData
{
    StaticArray<NodeID,2>                               nodeIndices;
    DynamicArray<DoFID>                                 dofMap;
    std::shared_ptr<TAnsatzSpace>                       pAnsatzSpace;
    std::shared_ptr<typename TAnsatzSpace::Derivative>  pAnsatzDerivatives;
}; // struct CellData


struct VertexData
{
    double position;
}; // struct VertexData


/** @brief 1D system test.
 *  @details This test models linear steady-state heat convection on a 1D bar
 *           with Neumann conditions on both boundaries. The heat flux on the
 *           left boundary is set to 1 while no heat flux is prescribed on the
 *           right boundary.
 *
 *           The mesh conforms to the boundaries and consists of 10 1D elements
 *           with linear basis functions.
 *
 *           Mesh:
 *           @code
 *                    +---+   +---+   +---+   +---+   +---+   +---+   +---+   +---+   +---+   +---+
 *           [q(0)=1] | 0 |-0-| 1 |-1-| 2 |-2-| 3 |-3-| 4 |-4-| 5 |-5-| 6 |-6-| 7 |-7-| 8 |-8-| 9 | [q(L)=0]
 *                    +---+   +---+   +---+   +---+   +---+   +---+   +---+   +---+   +---+   +---+
 *           @endcode
 */
CIE_TEST_CASE("1D", "[systemTests]")
{
    CIE_TEST_CASE_INIT("1D")
    constexpr unsigned Dimension = 1;
    constexpr Size nodeCount = 10;

    // Construct a linear 1D ansatz space
    using Basis = maths::Polynomial<Scalar>;
    using Ansatz = maths::AnsatzSpace<Basis,Dimension>;
    auto pAnsatzSpace = std::make_shared<Ansatz>(Ansatz::AnsatzSet {
        Basis({ 0.5,  0.5}),
        Basis({ 0.5, -0.5})
    });
    auto pAnsatzDerivatives = std::make_shared<Ansatz::Derivative>(pAnsatzSpace->makeDerivative());

    // Find ansatz functions that coincide on opposite boundaries.
    // In adjacent elements, these ansatz functions will have to map
    // to the same DoF in the assembled system.
    const StaticArray<Scalar,5> samples {-1.0, -0.5, 0.0, 0.5, 1.0};
    const auto ansatzMap = makeAnsatzMap(*pAnsatzSpace,
                                         samples,
                                         utils::Comparison<Scalar>(1e-8, 1e-6));

    // Construct mesh
    // The mesh is modeled by an adjacency graph whose vertices are cells,
    // and edges represent the boundaries between adjacent cells. The graph
    // is not directed, but the BoundaryID of the graph edge is always
    // defined on the edge's source cell.
    using Mesh = Graph<void,BoundaryID>;
    Mesh mesh;

    // Insert cells into the adjacency graph
    for (Size iCell : std::ranges::views::iota(0ul, nodeCount)) {
        // Insert the cell into the adjacency graph (mesh) as a vertex
        mesh.insert(Mesh::Vertex(
            Mesh::VertexID(iCell),
            {} // <== edges of the adjacency graph are added automatically during edge insertion
        ));
    }

    AttributeContainer<CellData<Ansatz>> cellAttributes(mesh.vertices().size());

    // Insert shared element boundaries into the adjacency graph as edges
    for (Size iBoundary : std::ranges::views::iota(0ul, nodeCount - 1)) {
        const Size iLeftCell = iBoundary;
        const Size iRightCell = iBoundary + 1;
        mesh.insert(Mesh::Edge(
            iBoundary,
            {iLeftCell, iRightCell},
            "+x" // <== all bundaries connect cells on the left to cells on the right
        ));
    }

    // Define and fill cell attributes
    std::vector<utils::SymmetricPair<Size>> dofPairs;
    for (Ref<const Mesh::Edge> rBoundary : mesh.edges()) {
        //Ref<CellData<Ansatz>> rSourceAttributes = cellAttributes.at<CellData<Ansatz>>(rBoundary.source());
        //Ref<CellData<Ansatz>> rTargetAttributes = cellAttributes.at<CellData<Ansatz>>(rBoundary.target());

        // Get coincident DoFs
        dofPairs.clear();
        ansatzMap.getPairs(rBoundary.data(), std::back_inserter(dofPairs));

        // Assign global DoF IDs to local DoFs. The following scenarios are possible:
        // - neither cells' DoFs have been assigned global DoF IDs yet. This case is the
        //   simplest, since they can be issued new ones and we only have to take special
        //   care of coincident DoFs.
        // - one of the cells has already been assigned global DoF IDs. This case is also
        //   straightforward, as coincident DoFs of the unassigned cell can be assigned
        //   global DoF IDs from the other cell.
        // - both cells have already been assigned global DoF IDs. This one's tricky,
        //   because the global IDs of coincident DoFs must match, which might require
        //   deleting existing DoFs if this condition is not already satisfied.
    }

    Assembler assembler;
    assembler.addGraph(mesh,
                       []([[maybe_unused]] Ref<const Mesh::Vertex> rVertex) {return 2ul;},
                       [&ansatzMap](Ref<const Mesh::Edge> rEdge, Assembler::DoFPairIterator it) {ansatzMap.getPairs(rEdge.data(), it);});

    // Create empty CSR matrix
    std::size_t rowCount, columnCount;
    DynamicArray<std::size_t> rowExtents, columnIndices;
    DynamicArray<double> nonzeros;
    assembler.makeCSRMatrix(rowCount, columnCount, rowExtents, columnIndices, nonzeros);

    {
        std::ofstream file("output.mm");
        utils::io::MatrixMarket::Output io(file);
        io(rowCount,
           columnCount,
           nonzeros.size(),
           rowExtents.data(),
           columnIndices.data(),
           nonzeros.data());
    }
}


} // namespace cie::fem
