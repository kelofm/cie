// --- Utility Includes ---
#include "packages/testing/inc/essentials.hpp"

// --- FEM Includes ---
#include "packages/graph/inc/Partition.hpp"
#include "packages/graph/inc/PartitionManager.hpp"
#include "packages/graph/inc/Graph.hpp"
#include "packages/graph/inc/BoundaryID.hpp"
#include "packages/maths/inc/LambdaExpression.hpp"
#include "packages/maths/inc/ScaleTranslateTransform.hpp"
#include "packages/numeric/inc/GaussLegendreQuadrature.hpp"
#include "packages/maths/inc/Polynomial.hpp"
#include "packages/maths/inc/AnsatzSpace.hpp"

// --- STL Includes ---
#include <ranges>
#include <iostream>


namespace cie::fem {


template <class TAnsatzSpace>
struct CellData
{
    StaticArray<Size,2>                                 nodeIndices;
    StaticArray<Size,2>                                 dofMap;
    std::shared_ptr<TAnsatzSpace>                       p_ansatzSpace;
    std::shared_ptr<typename TAnsatzSpace::Derivative>  p_ansatzDerivatives;
}; // struct CellData


struct VertexData
{
    double position;
};


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

    using Number = double;
    constexpr unsigned Dimension = 1;

    // Construct a linear 1D ansatz space
    using Basis = maths::Polynomial<Number>;
    using Ansatz = maths::AnsatzSpace<Basis,Dimension>;
    auto p_ansatzSpace = std::make_shared<Ansatz>(Ansatz::AnsatzSet {
        Basis({ 0.5,  0.5}),
        Basis({ 0.5, -0.5})
    });
    auto p_ansatzDerivatives = std::make_shared<Ansatz::Derivative>(p_ansatzSpace->makeDerivative());

    // Construct mesh
    // The mesh is modeled by a graph whose vertices are cells, and
    // edges represent the boundaries between adjacent cells. The graph
    // is not directed, but the BoundaryID of the graph edge is always
    // defined on the edge's source cell.
    using Mesh = Graph<void,BoundaryID>;
    Mesh mesh;
    for (Size i_element : std::ranges::views::iota(0, 9)) {
        const Size i_leftNode = i_element;
        const Size i_rightNode = i_element + 1;
        mesh.insert(Mesh::Edge(i_element, {i_leftNode, i_rightNode}));
    }

    // Set graph properties
    for (Mesh::Edge& r_boundary : mesh.edges()) {
        // All boundaries connect cells on the left to cells on the right,
        // so each edge's direction is in the positive x direction.
        r_boundary.data() = BoundaryID("+x");
    }

    // Construct cell attributes
    AttributeContainer<
        CellData<Ansatz>
    > cells;

    // Construct node attributes
    AttributeContainer<
        Number // <== position
    > nodes;


}


} // namespace cie::fem
