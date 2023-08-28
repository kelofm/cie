// --- FEM Includes ---
#include "packages/graph/inc/Partition.hpp"
#include "packages/graph/inc/PartitionManager.hpp"

// --- Utility Includes ---
#include "packages/testing/inc/essentials.hpp"
#include <packages/compile_time/packages/parameter_pack/inc/Sequence.hpp>


namespace cie::fem {


CIE_TEST_CASE("PartitionManager", "[graph]")
{
    CIE_TEST_CASE_INIT("PartitionManager")

    AttributeContainer<> empty;
    AttributeContainer<double> root;
    AttributeContainer<int, ParentIndex> mid;
    AttributeContainer<bool, ParentIndex> leaf;

    for (int i=0; i<10; ++i) {
        root.push_back(double(i) / 10.0);
    }

    for (std::size_t i=0; i<4; ++i) {
        mid.push_back(10 * i, 2 * i);
    }

    leaf.push_back(false, 3ul);
    leaf.push_back(true, 0ul);

    // root:   0.0   0.1   0.2   0.3   0.4   0.5   0.6   0.7   0.8   0.9
    // mid:    0           10          20          30
    // leaf:   true                                false

    Partition<decltype(root),decltype(empty)> rootPartition(std::move(root), decltype(empty)());
    Partition<decltype(mid),decltype(empty)> midPartition(std::move(mid), decltype(empty)());
    Partition<decltype(leaf),decltype(empty)> leafPartition(std::move(leaf), decltype(empty)());

    {
        const auto [d, i, b] = PartitionManager::collectAttributes<PartitionBase::Item::Vertex>(
            0,
            rootPartition,
            midPartition,
            leafPartition
        );
        CIE_TEST_CHECK(d == Approx(0.6));
        CIE_TEST_CHECK(i == 30);
        CIE_TEST_CHECK(b == false);
    }

    {
        const auto [d, i, b] = PartitionManager::collectAttributes<PartitionBase::Item::Vertex>(
            1,
            rootPartition,
            midPartition,
            leafPartition
        );
        CIE_TEST_CHECK(d == 0.0);
        CIE_TEST_CHECK(i == 0);
        CIE_TEST_CHECK(b == true);
    }
}


} // namespace cie::fem
