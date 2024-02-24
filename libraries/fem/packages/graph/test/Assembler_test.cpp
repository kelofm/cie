// --- Utility Includes ---
#include "packages/testing/inc/essentials.hpp"

// --- FEM Includes ---
#include "packages/graph/inc/Assembler.hpp"
#include "packages/maths/inc/Polynomial.hpp"
#include "packages/maths/inc/AnsatzSpace.hpp"
#include <iostream>


namespace cie::fem {


CIE_TEST_CASE("Assembler", "[graph]")
{
    CIE_TEST_CASE_INIT("Assembler")
    using VertexData = void;
    using EdgeData = BoundaryID;
    using Mesh = Graph<VertexData,EdgeData>;

    // Build mesh
    // + --- + --- + --- +
    // |  3  |  4  |  5  |
    // + --- + --- + --- +
    // |  0  |  1  |  2  |
    // + --- + --- + --- +
    Mesh mesh;
    std::size_t edgeID = 0ul;
    mesh.insert(Mesh::Edge(edgeID++, {0, 1}, "+x"));
    mesh.insert(Mesh::Edge(edgeID++, {1, 2}, "+x"));
    mesh.insert(Mesh::Edge(edgeID++, {3, 4}, "+x"));
    mesh.insert(Mesh::Edge(edgeID++, {4, 5}, "+x"));
    mesh.insert(Mesh::Edge(edgeID++, {0, 3}, "+y"));
    mesh.insert(Mesh::Edge(edgeID++, {1, 4}, "+y"));
    mesh.insert(Mesh::Edge(edgeID++, {2, 5}, "+y"));

    // Construct basis
    using Basis = maths::Polynomial<float>;
    using Ansatz = maths::AnsatzSpace<Basis,/*Dimension = */2>;
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
        ansatzMap.getPairs(rBoundary.data(), itOutput);
    };

    // Assign DoF indices
    Assembler assembler;
    assembler.addGraph(mesh, dofCounter, dofMatcher);

    // Check connectivities
    CIE_TEST_CHECK(assembler.dofCount() == 12);

    CIE_TEST_CHECK(assembler[0][1] == assembler[1][0]);
    CIE_TEST_CHECK(assembler[0][2] == assembler[3][0]);
    CIE_TEST_CHECK(assembler[0][3] == assembler[1][2]);
    CIE_TEST_CHECK(assembler[0][3] == assembler[3][1]);
    CIE_TEST_CHECK(assembler[0][3] == assembler[4][0]);

    CIE_TEST_CHECK(assembler[1][1] == assembler[2][0]);
    CIE_TEST_CHECK(assembler[1][2] == assembler[3][1]);
    CIE_TEST_CHECK(assembler[1][2] == assembler[4][0]);
    CIE_TEST_CHECK(assembler[1][3] == assembler[2][2]);
    CIE_TEST_CHECK(assembler[1][3] == assembler[4][1]);
    CIE_TEST_CHECK(assembler[1][3] == assembler[5][0]);

    CIE_TEST_CHECK(assembler[2][2] == assembler[4][1]);
    CIE_TEST_CHECK(assembler[2][2] == assembler[5][0]);
    CIE_TEST_CHECK(assembler[2][3] == assembler[5][1]);
} // CIE_TEST_CASE "Assembler"


} // namespace cie::fem
