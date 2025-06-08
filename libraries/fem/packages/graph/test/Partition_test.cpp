// --- FEM Includes ---
#include "packages/graph/inc/Partition.hpp"
#include "packages/graph/inc/PartitionVisitor.hpp"
#include "packages/maths/inc/Expression.hpp"
#include "packages/maths/inc/LambdaExpression.hpp"
#include "packages/maths/inc/ScaleTranslateTransform.hpp"
#include "packages/numeric/inc/GaussLegendreQuadrature.hpp"
#include "packages/maths/inc/Polynomial.hpp"
#include "packages/maths/inc/AnsatzSpace.hpp"

// --- Utility Includes ---
#include "packages/numeric/inc/Quadrature.hpp"
#include "packages/testing/inc/essentials.hpp"
#include "packages/stl_extension/inc/StaticArray.hpp"
#include "packages/matrix/inc/DynamicEigenMatrix.hpp"
#include "packages/ranges/inc/TransformIterator.hpp"

// --- STL Includes ---
#include <ranges>
#include <iostream>


namespace cie::fem {


CIE_TEST_CASE("PartitionVisitor", "[graph]")
{
    CIE_TEST_CASE_INIT("PartitionVisitor")

    AttributeContainer<double> root;
    AttributeContainer<int, unsigned, ParentIndex> mid;
    AttributeContainer<char, ParentIndex> leaf;

    for (int i=0; i<10; ++i) {
        root.push_back(double(i) / 10.0);
    }
    CIE_TEST_CHECK(root.size() == 10);

    for (std::size_t i=0; i<4; ++i) {
        mid.push_back(10 * i, 5u * i, 2ul * i);
    }
    CIE_TEST_CHECK(mid.size() == 4);

    leaf.push_back(false, 3ul);
    leaf.push_back(true, 0ul);
    CIE_TEST_CHECK(leaf.size() == 2);

    // root:   0.0   0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9
    // mid:    0           10          20          30
    // leaf:   true                                false

    using LeafPartition = PartitionTree<Partition<decltype(leaf)>>;
    using MidPartition = PartitionTree<Partition<decltype(mid)>,LeafPartition>;
    using RootPartition = PartitionTree<Partition<decltype(root)>,MidPartition>;

    const auto pRootPartition = std::make_shared<RootPartition>("root", 0, std::move(root));
    {
        auto pMidPartition = std::make_unique<MidPartition>("mid", 1, std::move(mid));
        auto pLeafPartition = std::make_unique<LeafPartition>("leaf", 2, std::move(leaf));

        pRootPartition->setChild<0>(std::move(pMidPartition));
        pRootPartition->setChild<0,0>(std::move(pLeafPartition));
    }

    // Collect attributes from a set of partition references
    {
        const auto& [d, i, u, b] = PartitionVisitor::collectCellAttributes(
            0,
            *pRootPartition,
            pRootPartition->child<0>(),
            pRootPartition->child<0,0>());
        CIE_TEST_CHECK(d == Approx(0.6));
        CIE_TEST_CHECK(i == 30);
        CIE_TEST_CHECK(u == 15u);
        CIE_TEST_CHECK(b == false);
    }

    {
        const auto [d, i, u, b] = PartitionVisitor::collectCellAttributes(
            1,
            *pRootPartition,
            pRootPartition->child<0>(),
            pRootPartition->child<0,0>());
        CIE_TEST_CHECK(d == 0.0);
        CIE_TEST_CHECK(i == 0);
        CIE_TEST_CHECK(u == 0u);
        CIE_TEST_CHECK(b == true);
    }

    // Collect attributes from a set of partition pointers
    {
        const auto [d, i, u, b] = PartitionVisitor::collectCellAttributes(
            0,
            &*pRootPartition,
            &pRootPartition->child<0>(),
            &pRootPartition->child<0,0>());
        CIE_TEST_CHECK(d == Approx(0.6));
        CIE_TEST_CHECK(i == 30);
        CIE_TEST_CHECK(u == 15u);
        CIE_TEST_CHECK(b == false);
    }

    {
        const auto [d, i, u, b] = PartitionVisitor::collectCellAttributes(
            1,
            &*pRootPartition,
            &pRootPartition->child<0>(),
            &pRootPartition->child<0,0>());
        CIE_TEST_CHECK(d == 0.0);
        CIE_TEST_CHECK(i == 0);
        CIE_TEST_CHECK(u == 0u);
        CIE_TEST_CHECK(b == true);
    }

    // Collect attributes from a tuple of attribute pointers.
    {
        const auto [d, i, u, b] = PartitionVisitor::collectCellAttributes(
            0,
            std::make_tuple(&*pRootPartition,
                            &pRootPartition->child<0>(),
                            &pRootPartition->child<0,0>()));
        CIE_TEST_CHECK(d == Approx(0.6));
        CIE_TEST_CHECK(i == 30);
        CIE_TEST_CHECK(u == 15u);
        CIE_TEST_CHECK(b == false);
    }

    {
        const auto [d, i, u, b] = PartitionVisitor::collectCellAttributes(
            1,
            std::make_tuple(&*pRootPartition,
                            &pRootPartition->child<0>(),
                            &pRootPartition->child<0,0>()));
        CIE_TEST_CHECK(d == 0.0);
        CIE_TEST_CHECK(i == 0);
        CIE_TEST_CHECK(u == 0u);
        CIE_TEST_CHECK(b == true);
    }

    // Collect attributes from a tuple of attribute pointers
    // obtained from flattening a partition tree.
    {
        const auto [d, i, u, b] = PartitionVisitor::collectCellAttributes(
            0,
            makePartitionHandle<0,0>(*pRootPartition).flatten());
        CIE_TEST_CHECK(d == Approx(0.6));
        CIE_TEST_CHECK(i == 30);
        CIE_TEST_CHECK(u == 15u);
        CIE_TEST_CHECK(b == false);
    }

    {
        const auto [d, i, u, b] = PartitionVisitor::collectCellAttributes(
            1,
            makePartitionHandle<0,0>(*pRootPartition).flatten());
        CIE_TEST_CHECK(d == 0.0);
        CIE_TEST_CHECK(i == 0);
        CIE_TEST_CHECK(u == 0u);
        CIE_TEST_CHECK(b == true);
    }

    // Collect mutable attributes.
    {
        auto partitionHandle = makePartitionHandle<0,0>(*pRootPartition);
        const auto [d, i, u, b] = PartitionVisitor::collectCellAttributes(
            0,
            partitionHandle.flatten());
        d = 2.0;
        CIE_TEST_CHECK(
            std::get<0>(PartitionVisitor::collectCellAttributes(0, partitionHandle.flatten()))
            ==
            Approx(2.0)
        );
    }

    // Collect attributes from a partial tree.
    {
        const auto [i, u, b] = PartitionVisitor::collectCellAttributes(
            0,
            makePartitionHandle<0>(pRootPartition->child<0>()).flatten());
        CIE_TEST_CHECK(i == 30);
        CIE_TEST_CHECK(u == 15u);
        CIE_TEST_CHECK(b == false);
    }

    {
        const auto [i, u, b] = PartitionVisitor::collectCellAttributes(
            1,
            makePartitionHandle<0>(pRootPartition->child<0>()).flatten());
        CIE_TEST_CHECK(i == 0);
        CIE_TEST_CHECK(u == 0u);
        CIE_TEST_CHECK(b == true);
    }

    // Collect attributes from an immutable partial tree.
    {
        const auto partitionHandle = makePartitionHandle<0>(pRootPartition->child<0>());
        const auto& [i, u, b] = PartitionVisitor::collectCellAttributes(
            0,
            partitionHandle.flatten());
        CIE_TEST_CHECK(i == 30);
        CIE_TEST_CHECK(u == 15u);
        CIE_TEST_CHECK(b == false);
    }

    {
        const auto partitionHandle = makePartitionHandle<0>(pRootPartition->child<0>());
        const auto [i, u, b] = PartitionVisitor::collectCellAttributes(
            1,
            partitionHandle.flatten());
        CIE_TEST_CHECK(i == 0);
        CIE_TEST_CHECK(u == 0u);
        CIE_TEST_CHECK(b == true);
    }

    // Collect and mutate attributes from a mutable partial tree.
    {
        auto [i, u, b] = PartitionVisitor::collectCellAttributes(
            0,
            makePartitionHandle<0>(pRootPartition->child<0>()).flatten());

        CIE_TEST_CHECK(i == 30);
        CIE_TEST_CHECK(u == 15u);
        CIE_TEST_CHECK(b == false);

        i *= 2;
        u *= 2;
        b = !b;

        CIE_TEST_CHECK(i == 60);
        CIE_TEST_CHECK(u == 30u);
        CIE_TEST_CHECK(b == true);

        const auto& [ii, uu, bb] = PartitionVisitor::collectCellAttributes(
            0,
            makePartitionHandle<0>(pRootPartition->child<0>()).flatten());

        CIE_TEST_CHECK(i == ii);
        CIE_TEST_CHECK(u == uu);
        CIE_TEST_CHECK(b == bb);

        CIE_TEST_CHECK(ii == 60);
        CIE_TEST_CHECK(uu == 30u);
        CIE_TEST_CHECK(bb == true);
    }

    {
        auto [i, u, b] = PartitionVisitor::collectCellAttributes(
            1,
            makePartitionHandle<0>(pRootPartition->child<0>()).flatten());

        CIE_TEST_CHECK(i == 0);
        CIE_TEST_CHECK(u == 0u);
        CIE_TEST_CHECK(b == true);

        i += 2;
        u += 3;
        b = !b;

        CIE_TEST_CHECK(i == 2);
        CIE_TEST_CHECK(u == 3u);
        CIE_TEST_CHECK(b == false);

        const auto& [ii, uu, bb] = PartitionVisitor::collectCellAttributes(
            1,
            makePartitionHandle<0>(pRootPartition->child<0>()).flatten());

        CIE_TEST_CHECK(i == ii);
        CIE_TEST_CHECK(u == uu);
        CIE_TEST_CHECK(b == bb);

        CIE_TEST_CHECK(ii == 2);
        CIE_TEST_CHECK(uu == 3u);
        CIE_TEST_CHECK(bb == false);
    }
}


} // namespace cie::fem
