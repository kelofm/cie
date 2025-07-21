// --- Utility Includes ---
#include "packages/testing/inc/essentials.hpp"
#include "packages/io/inc/MatrixMarket.hpp"

// --- FEM Includes ---
#include "packages/graph/inc/Assembler.hpp"
#include "packages/maths/inc/Polynomial.hpp"
#include "packages/maths/inc/AnsatzSpace.hpp"
#include "packages/graph/inc/OrientedBoundary.hpp"
#include "packages/graph/inc/connectivity.hpp"
#include <packages/io/inc/MatrixMarket.hpp>
#include <packages/macros/inc/exceptions.hpp>


namespace cie::fem {


CIE_TEST_CASE("Assembler", "[graph]")
{
    CIE_TEST_CASE_INIT("Assembler")
    using Boundary = OrientedBoundary<2>;
    using VertexData = void;
    using EdgeData = std::pair<Boundary,Boundary>;
    using Mesh = Graph<VertexData,EdgeData>;

    // Build mesh
    /** Build a mesh that consists of 6 linear quads with a linear set
     *  of Lagrangian basis functions. The local axes of each quad are
     *  flipped such that functions that don't vanish on boundaries
     *  coincide with the corresponding functions of neighboring cells.
     *
     *           0 ----------- 1                           1 ----------- 0                           0 ----------- 1
     *           |             |                           |             |                           |             |
     *           |             |                           |             |                           |             |
     *           |      3 > x  | "+x-y:+x"       "-x-y:+x" |  x < 4      | "-x-y:-x"       "+x-y:-x" |      5 > x  |
     *           |      v      |                           |      v      |                           |      v      |
     *           |      y      |                           |      y      |                           |      y      |
     *           2 ----------- 3                           3 ----------- 2                           2 ----------- 3
     *              "+x-y:+y"                                 "-x-y:+y"                                 "-x-y:+y"
     *
     *              "+x+y:+y"                                 "-x+y:+y"                                 "-x+y:+y"
     *           2 ----------- 3                           3 ----------- 2                           2 ----------- 3
     *           |      y      |                           |      y      |                           |      y      |
     *           |      ^      |                           |      ^      |                           |      ^      |
     *           |      0 > x  | "+x+y:+x"       "-x+y:+x" |  x < 1      | "-x+y:-x"       "+x+y:-x" |      2 > x  |
     *           |             |                           |             |                           |             |
     *           |             |                           |             |                           |             |
     *           0 ----------- 1                           1 ----------- 0                           0 ----------- 1
     */
    Mesh mesh;
    std::size_t edgeID = 0ul;
    mesh.insert(Mesh::Edge(edgeID++,
                           {0, 1},
                           std::make_pair(Boundary("+x+y", "+x"),
                                          Boundary("-x+y", "+x"))));
    mesh.insert(Mesh::Edge(edgeID++,
                           {1, 2},
                           std::make_pair(Boundary("-x+y", "-x"),
                                          Boundary("+x+y", "-x"))));
    mesh.insert(Mesh::Edge(edgeID++,
                           {3, 4},
                           std::make_pair(Boundary("+x-y", "+x"),
                                          Boundary("-x-y", "+x"))));
    mesh.insert(Mesh::Edge(edgeID++,
                           {4, 5},
                           std::make_pair(Boundary("-x-y", "-x"),
                                          Boundary("+x-y", "-x"))));
    mesh.insert(Mesh::Edge(edgeID++,
                           {0, 3},
                           std::make_pair(Boundary("+x+y", "+y"),
                                          Boundary("+x-y", "+y"))));
    mesh.insert(Mesh::Edge(edgeID++,
                           {1, 4},
                           std::make_pair(Boundary("-x+y", "+y"),
                                          Boundary("-x-y", "+y"))));
    mesh.insert(Mesh::Edge(edgeID++,
                           {2, 5},
                           std::make_pair(Boundary("+x+y", "+y"),
                                          Boundary("+x-y", "+y"))));

    // Construct basis
    using Basis = maths::Polynomial<float>;
    using Ansatz = maths::AnsatzSpace<Basis,/*Dimension=*/2>;
    const auto pAnsatzSpace = std::make_shared<Ansatz>(Ansatz::AnsatzSet {
        Basis({ 0.5, -0.5}),
        Basis({ 0.5,  0.5})
    });

    const auto ansatzMap = makeAnsatzMap(*pAnsatzSpace,
                                         DynamicArray<float> {-1.0f, -0.5f, 0.0f, 0.5f, 1.0f},
                                         utils::Comparison<float>(1e-4, 1e-3));
    const auto dofCounter = [&pAnsatzSpace]([[maybe_unused]] const auto& _) -> std::size_t {
        return pAnsatzSpace->size();
    };
    const auto dofMatcher = [&ansatzMap](Ref<const Mesh::Edge> rBoundary,
                                         Assembler::DoFPairIterator itOutput) -> void {
        CIE_BEGIN_EXCEPTION_TRACING
        ansatzMap.getPairs(rBoundary.data().first,
                           rBoundary.data().second,
                           itOutput);
        CIE_END_EXCEPTION_TRACING
    };

    // Assign DoF indices
    Assembler assembler;
    assembler.addGraph(mesh, dofCounter, dofMatcher);

    // Check connectivities
    CIE_TEST_CHECK(assembler.dofCount() == 12);

    CIE_TEST_CHECK(assembler[0][1] == assembler[1][1]);
    CIE_TEST_CHECK(assembler[0][2] == assembler[3][2]);
    CIE_TEST_CHECK(assembler[0][3] == assembler[1][3]);
    CIE_TEST_CHECK(assembler[0][3] == assembler[3][3]);
    CIE_TEST_CHECK(assembler[0][3] == assembler[4][3]);

    CIE_TEST_CHECK(assembler[1][0] == assembler[2][0]);
    CIE_TEST_CHECK(assembler[1][2] == assembler[2][2]);
    CIE_TEST_CHECK(assembler[1][2] == assembler[4][2]);
    CIE_TEST_CHECK(assembler[1][2] == assembler[5][2]);

    CIE_TEST_CHECK(assembler[2][2] == assembler[4][2]);
    CIE_TEST_CHECK(assembler[2][2] == assembler[5][2]);
    CIE_TEST_CHECK(assembler[2][3] == assembler[5][3]);

    CIE_TEST_CHECK(assembler[3][1] == assembler[4][1]);
    CIE_TEST_CHECK(assembler[4][0] == assembler[5][0]);

    {
        int rowCount, columnCount;
        DynamicArray<int> rowExtents, columnIndices;
        DynamicArray<float> entries;
        assembler.makeCSRMatrix(rowCount,
                                columnCount,
                                rowExtents,
                                columnIndices,
                                entries);
        std::ofstream file("assembler_test_2d.mm");
        cie::io::MatrixMarket::Output output(file);
        CIE_TEST_CHECK_NOTHROW(output(rowCount,
                                      columnCount,
                                      entries.size(),
                                      rowExtents.data(),
                                      columnIndices.data(),
                                      entries.data()));
    }
} // CIE_TEST_CASE "Assembler"


} // namespace cie::fem
