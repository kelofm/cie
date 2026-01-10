#define _USE_MATH_DEFINES

// --- Utility Includes ---
#include "packages/testing/inc/essentials.hpp"
#include "packages/concurrency/inc/sycl.hpp"

// --- Internal Includes ---
#include "packages/partitioning/inc/AABBoxNode.hpp"

// --- STL Inlcudes ---
#include <cmath>


namespace cie::geo {


class TestAABBoxNodeObjectType final
    : public BoxBoundable<2,Double>,
      public AABBox<2,Double>
{
public:
    using AABBox<2,Double>::Point;

    TestAABBoxNodeObjectType(Ref<const Point> rBase,
                             Ref<const Point> rLengths)
        : BoxBoundable<2,Double>(),
          AABBox<2,Double>(rBase, rLengths) {}
private:
    void computeBoundingBoxImpl(TestAABBoxNodeObjectType::BoundingBox& rBox) noexcept override
    {rBox = *this;}
};


CIE_TEST_CASE("AABBoxNode", "[partitioning]")
{
    CIE_TEST_CASE_INIT("AABBoxNode")

    using TCoordinate       = TestAABBoxNodeObjectType::Coordinate;
    using Object            = TestAABBoxNodeObjectType;
    using Node              = AABBoxNode<Object>;

    // Generate objects
    std::vector<Object> objects;

    constexpr Size cellsPerDirection = 5;

    for (Size i=0; i<cellsPerDirection; ++i) {
        for (Size j=0; j<cellsPerDirection; ++j) {
            if (i > j) {
                objects.emplace_back(
                    Object::Point {double(i)/cellsPerDirection, double(j)/cellsPerDirection},
                    Object::Point {1.0 / cellsPerDirection, 1.0 / cellsPerDirection});
            } else {
                objects.emplace_back(
                    Object::Point {std::numeric_limits<TCoordinate>::max(), std::numeric_limits<TCoordinate>::max()},
                    Object::Point {0.0, 0.0}
                );
            }
        }
    }

    // Generate root node
    constexpr double delta = 0.0;
    Node root(
        Object::Point {-delta, -delta},
        Object::Point {1.0+2*delta, 1.0+2*delta},
        nullptr);

    for (auto& rObject : objects)
        CIE_TEST_CHECK_NOTHROW(root.insert(&rObject));

    CIE_TEST_CHECK(root.contained().size() == cellsPerDirection * (cellsPerDirection - 1) / 2);

    // Partition
    constexpr Size maxObjects = 3;
    constexpr Size maxLevel   = 5;
    bool partitionSuccess = false;

    CIE_TEST_CHECK_NOTHROW(partitionSuccess = root.partition(maxObjects, maxLevel));
    CIE_TEST_CHECK(partitionSuccess);

    CIE_TEST_CHECK_NOTHROW(root.shrink());

    // Check number of objects and maximum levels
    {
        auto nodeVisitFunction = [maxLevel=maxLevel,maxObjects=maxObjects](Node* p_node) -> bool {
            CIE_TEST_CHECK(p_node->level() <= maxLevel);
            if (p_node->isLeaf())
                CIE_TEST_CHECK(p_node->contained().size() <= maxObjects);
            return true;
        };
        CIE_TEST_CHECK_NOTHROW(root.visit(nodeVisitFunction));
    }
}


CIE_TEST_CASE("FlatAABBoxNode", "[partitioning]")
{
    CIE_TEST_CASE_INIT("FlatAABBoxNode")

    const Size Dimension    = TestAABBoxNodeObjectType::Dimension;
    using TCoordinate       = TestAABBoxNodeObjectType::Coordinate;
    using Object            = TestAABBoxNodeObjectType;
    using Node              = AABBoxNode<Object>;

    // Generate objects
    std::vector<Object> objects;

    constexpr Size cellsPerDirection = 5;

    for (Size i=0; i<cellsPerDirection; ++i) {
        for (Size j=0; j<cellsPerDirection; ++j) {
            if (i > j) {
                objects.emplace_back(
                    Object::Point {double(i)/cellsPerDirection, double(j)/cellsPerDirection},
                    Object::Point {1.0 / cellsPerDirection, 1.0 / cellsPerDirection});
            } else {
                objects.emplace_back(
                    Object::Point {std::numeric_limits<TCoordinate>::max(), std::numeric_limits<TCoordinate>::max()},
                    Object::Point {0.0, 0.0}
                );
            }
        }
    }

    // Generate root node
    double delta = 0.0;
    Node root(
        Object::Point {-delta, -delta},
        Object::Point {1.0+2*delta, 1.0+2*delta},
        nullptr);

    for (auto& rObject : objects)
        root.insert(&rObject);

    // Partition
    constexpr Size maxObjects = 3;
    constexpr Size maxLevel   = 5;
    root.partition(maxObjects, maxLevel);
    root.shrink();

    const auto flatTree = FlatAABBoxTree<
        TCoordinate,
        Dimension,
        unsigned,
        std::allocator<std::byte>
    >::flatten(
        root,
        [&objects] (Ref<const Object> rObject) -> unsigned {
            return std::distance(static_cast<Ptr<const Object>>(objects.data()), &rObject);
        },
        std::allocator<std::byte>());

    for (unsigned i=0; i<cellsPerDirection; ++i) {
        for (unsigned j=0; j<cellsPerDirection; ++j) {
            const Object::Point center {
                (double(i) + 0.5) / cellsPerDirection,
                (double(j) + 0.5) / cellsPerDirection
            };

            CIE_TEST_CHECK_NOTHROW(flatTree.makeView().find(
                std::span<const TCoordinate,Dimension>(center.data(), Dimension),
                std::span<const Object> {objects.data(), objects.size()}
            ));

            const std::size_t iMaybeObject = flatTree.makeView().find(
                std::span<const TCoordinate,Dimension>(center.data(), Dimension),
                std::span<const Object> {objects.data(), objects.size()}
            );

            if (i > j) {
                CIE_TEST_REQUIRE(iMaybeObject != objects.size());
                CIE_TEST_CHECK(iMaybeObject == i * cellsPerDirection + j);
            } else {
                CIE_TEST_CHECK(iMaybeObject == objects.size());
            }
        }
    }
}


CIE_TEST_CASE("AABBoxNode with point objects", "[partitioning]") {
    CIE_TEST_CASE_INIT("AABBoxNode with point objects")

    const Size Dimension    = TestAABBoxNodeObjectType::Dimension;
    using TCoordinate       = TestAABBoxNodeObjectType::Coordinate;
    using TPoint            = typename Traits<Dimension,TCoordinate>::Point;
    using Object            = TPoint;
    using Node              = AABBoxNode<Object>;

    // Create points on a circle
    const Size numberOfObjects = 100;
    const double dt = 2.0 * M_PI / (numberOfObjects-1);

    std::vector<Object> objects;
    for (double t=0.0; t<2.0*M_PI; t+=dt)
        objects.emplace_back(Object {
            0.5 * std::cos(t) + 0.5,
            0.5 * std::sin(t) + 0.5 + 1e-10});

    // Construct root and add points to it
    constexpr double delta = 0.0;
    CIE_TEST_REQUIRE_NOTHROW(Node(
        Object {-delta, -delta},
        Object {1.0+2*delta, 1.0+2*delta},
        nullptr
    ));
    Node root(
        Object {-delta, -delta},
        Object {1.0+2*delta, 1.0+2*delta},
        nullptr
    );

    for (auto& rObject : objects)
        CIE_TEST_CHECK_NOTHROW(root.insert(&rObject));

    CIE_TEST_CHECK(root.contained().size() == numberOfObjects);
    CIE_TEST_CHECK(root.intersected().size() == 0);

    // Check partitioning
    constexpr Size maxObjects = 5;
    constexpr Size maxLevel   = 3;
    bool partitionSuccess = false;

    CIE_TEST_CHECK_NOTHROW(partitionSuccess = root.partition(maxObjects,maxLevel));
    CIE_TEST_CHECK(partitionSuccess);

    // Check number of objects and maximum level
    Size objectCounter = 0;

    auto nodeVisitFunctor = [&objectCounter, maxLevel=maxLevel, maxObjects=maxObjects](Node* p_node) -> bool {
        CIE_TEST_CHECK(p_node->intersected().size() == 0);
        CIE_TEST_CHECK(p_node->level() <= maxLevel);
        if (p_node->isLeaf()) {
            Size numberOfcontained = p_node->contained().size();
            CIE_TEST_CHECK(numberOfcontained <= maxObjects);
            objectCounter += numberOfcontained;
        }
        return true;
    };

    CIE_TEST_CHECK_NOTHROW(root.visit(nodeVisitFunctor));
    CIE_TEST_CHECK(objectCounter == numberOfObjects);
}


#ifdef CIE_ENABLE_SYCL
CIE_TEST_CASE("AABBoxNode SYCL", "[partitioning,sycl]")
{
    CIE_TEST_CASE_INIT("AABBoxNode SYCL")

    const Size Dimension = TestAABBoxNodeObjectType::Dimension;
    using TCoordinate    = TestAABBoxNodeObjectType::Coordinate;

    using Object    = stack::Box<Dimension,TCoordinate>;
    using Node      = AABBoxNode<Object>;

    constexpr Size cellsPerDirection = 5;
    constexpr std::size_t samplesPerDirection = 1e3;

    // SYCL setup.
    CIE_TEST_REQUIRE_NOTHROW(sycl::device(sycl::default_selector_v));
    sycl::device device(sycl::default_selector_v);
    sycl::queue queue(device);

    // Construct SYCL allocators.
    sycl::usm_allocator<Object::Point,sycl::usm::alloc::shared> pointAllocator(queue);
    sycl::usm_allocator<Object,sycl::usm::alloc::shared> objectAllocator(queue);
    sycl::usm_allocator<std::size_t,sycl::usm::alloc::shared> unsignedAllocator(queue);

    // Generate objects
    Ptr<Object> pObjectBegin = objectAllocator.allocate(cellsPerDirection * cellsPerDirection);
    for (Size i=0; i<cellsPerDirection; ++i) {
        for (Size j=0; j<cellsPerDirection; ++j) {
            if (i > j) {
                pObjectBegin[i * cellsPerDirection + j] = Object(
                    Object::Point {double(i)/cellsPerDirection, double(j)/cellsPerDirection},
                    Object::Point {1.0 / cellsPerDirection, 1.0 / cellsPerDirection});
            } else {
                pObjectBegin[i * cellsPerDirection + j] = Object(
                    Object::Point {std::numeric_limits<TCoordinate>::max(), std::numeric_limits<TCoordinate>::max()},
                    Object::Point {0.0, 0.0}
                );
            }
        }
    }

    // Generate root node.
    constexpr double delta = 0.0;
    Node root(
        Object::Point {-delta,          -delta},
        Object::Point {1.0 + 2 * delta, 1.0 + 2 * delta},
        nullptr);

    for (unsigned iObject=0u; iObject<cellsPerDirection*cellsPerDirection; ++iObject) {
        CIE_TEST_REQUIRE_NOTHROW(root.insert(pObjectBegin + iObject));
    }

    // Partition.
    constexpr Size maxObjects = 3;
    constexpr Size maxLevel   = 5;
    root.partition(maxObjects, maxLevel);
    root.shrink();

    // Flatten the tree into a shared memory region.
    const auto flatTree = FlatAABBoxTree<
        TCoordinate,
        Dimension,
        unsigned,
        sycl::usm_allocator<std::byte,sycl::usm::alloc::shared>
    >::flatten(
        root,
        [pObjectBegin] (Ref<const Object> rObject) -> unsigned {
            return std::distance(static_cast<Ptr<const Object>>(pObjectBegin), &rObject);
        },
        sycl::usm_allocator<std::byte,sycl::usm::alloc::shared>(queue)
    );

    // Construct sample points.
    Ptr<Object::Point> pSampleBegin = pointAllocator.allocate(samplesPerDirection * samplesPerDirection);
    queue.parallel_for(sycl::range(samplesPerDirection, samplesPerDirection), [pSampleBegin] (sycl::id<2> index) {
        Ref<Object::Point> rPoint = pSampleBegin[index.get(0) * samplesPerDirection + index.get(1)];
        rPoint[0] = (double(index.get(0)) + 0.5) / samplesPerDirection;
        rPoint[1] = (double(index.get(1)) + 0.5) / samplesPerDirection;
    }).wait();

    // Find which box each sample point is located in.
    Ptr<std::size_t> pMaybeObjectIndexBegin = unsignedAllocator.allocate(samplesPerDirection * samplesPerDirection);
    queue.parallel_for(
        sycl::range(samplesPerDirection * samplesPerDirection),
        [pMaybeObjectIndexBegin, pSampleBegin, pObjectBegin, treeView = flatTree.makeView()] (sycl::id<1> index) {
            pMaybeObjectIndexBegin[index.get(0)] = treeView.find(
                std::span<const TCoordinate,Dimension>(pSampleBegin[index.get(0)].data(), Dimension),
                std::span<const Object>(pObjectBegin, cellsPerDirection * cellsPerDirection));
    }).wait();

    // Check results.
    for (unsigned i=0; i<cellsPerDirection; ++i) {
        for (unsigned j=0; j<cellsPerDirection; ++j) {
            const Object::Point center {
                (double(i) + 0.5) / cellsPerDirection,
                (double(j) + 0.5) / cellsPerDirection
            };

            CIE_TEST_CHECK_NOTHROW(flatTree.makeView().find(
                std::span<const TCoordinate,Dimension>(center.data(), Dimension),
                std::span<const Object> {pObjectBegin, cellsPerDirection * cellsPerDirection}
            ));

            const std::size_t iMaybeObject = flatTree.makeView().find(
                std::span<const TCoordinate,Dimension>(center.data(), Dimension),
                std::span<const Object> {pObjectBegin, cellsPerDirection * cellsPerDirection}
            );

            if (i > j) {
                CIE_TEST_REQUIRE(iMaybeObject != cellsPerDirection * cellsPerDirection);
                CIE_TEST_CHECK(iMaybeObject == i * cellsPerDirection + j);
            } else {
                CIE_TEST_CHECK(iMaybeObject == cellsPerDirection * cellsPerDirection);
            }
        }
    }

    // Deallocate shared and device memory.
    unsignedAllocator.deallocate(pMaybeObjectIndexBegin, samplesPerDirection * samplesPerDirection);
    objectAllocator.deallocate(pObjectBegin, cellsPerDirection * cellsPerDirection);
    pointAllocator.deallocate(pSampleBegin, samplesPerDirection * samplesPerDirection);
}
#endif


} // namespace cie::geo
