// --- External Includes ---
#include "benchmark/benchmark.h"

// --- Utility Includes ---
#include "packages/maths/inc/Comparison.hpp"

// --- Geo Includes ---
#include "packages/trees/inc/ContiguousSpaceTree.hpp"
#include "packages/trees/inc/SpaceTreeNode.hpp"
#include <packages/primitives/inc/Cube.hpp>
#include <packages/primitives/inc/Box.hpp>
#include <packages/trees/inc/CartesianGridSampler.hpp>
#include <packages/trees/inc/MidPointSplitPolicy.hpp>
#include "packages/trees/inc/Cell.hpp"

// --- FEM Includes ---
#include "packages/maths/inc/OuterProduct.hpp"

// --- STL Includes ---
#include <array>
#include <atomic>
#include <limits>
#include <numbers>


// Deepest level of the tree
constexpr unsigned MAX_DEPTH = 8;

// Number of sample points in a cell per dimension
constexpr unsigned CELL_RESOLUTION = 5;


template <cie::concepts::Cube TGeometry>
TGeometry makeRootGeometry()
{
    typename TGeometry::Point base;
    std::fill(base.begin(),
              base.end(),
              typename TGeometry::Coordinate(-2));
    return TGeometry(base, typename TGeometry::Point::value_type(4));
}


template <cie::concepts::Box TGeometry>
TGeometry makeRootGeometry()
{
    typename TGeometry::Point base, lengths;
    std::fill(base.begin(),
              base.end(),
              typename TGeometry::Coordinate(-2));
    std::fill(lengths.begin(),
              lengths.end(),
              typename TGeometry::Coordinate(4));
    return TGeometry(base, lengths);
}


template <class TValue>
bool isInside(const TValue* it_begin,
              const TValue* it_end)
{
    TValue level = 1;
    for (; it_begin!=it_end; ++it_begin) {
        level *= (*it_begin) * (*it_begin);
    }
    return 0 < level;
}


template <class TGeometry>
void contiguousTree(benchmark::State& r_state)
{
    constexpr auto Dimension = TGeometry::Dimension;
    using Value = typename TGeometry::Coordinate;
    using Tree = cie::geo::ContiguousSpaceTree<TGeometry,std::size_t>;

    const cie::utils::Comparison<Value> comparison(
        std::numeric_limits<Value>::min(),
        5e-1
    );
    const Value referenceVolume =
        Dimension == 1 ? 2 :
        Dimension == 2 ? std::numbers::pi :
        Dimension == 3 ? 4.0 / 3.0 * std::numbers::pi :
        std::numeric_limits<Value>::max();

    for ([[maybe_unused]] auto dummy : r_state) {
        Tree tree(makeRootGeometry<TGeometry>());
        std::atomic<Value> volume = 0;

        const auto visitor = [&tree, &volume](const typename Tree::Node& r_node, unsigned depth) -> bool {
            if (MAX_DEPTH < depth) {return false;}

            // Query node geometry
            typename Tree::Point base, lengths;
            tree.getNodeGeometry(r_node, base.data(), lengths.data());
            if constexpr (cie::concepts::Cube<TGeometry>) {
                std::fill(lengths.begin() + 1,
                          lengths.end(),
                          lengths.front());
            }

            std::transform(lengths.begin(),
                           lengths.end(),
                           lengths.begin(),
                           [](Value length){return length / (CELL_RESOLUTION - 1);});

            // Loop over sample points
            std::array<std::size_t,Dimension> state;
            std::fill(state.begin(), state.end(), 0ul);
            unsigned char result = 0b1; // first bit indicates whether the result is uninitialized,
                                        // second bit indicates whether the sample is inside the target
                                        // third bit indicates whether the cell is split
            do {
                typename Tree::Point sample;
                for (unsigned i_dim=0; i_dim<Dimension; ++i_dim) {
                    sample[i_dim] = base[i_dim] + state[i_dim] * lengths[i_dim];
                }
                const bool inside = isInside(sample.data(), sample.data() + Dimension);
                if (result & 0b1) {
                    result = char(inside) << 1;
                } else if (inside != bool(result & 0b10)) {
                    result = 0b100;
                    break;
                }
            } while (cie::fem::maths::OuterProduct<Dimension>::next(CELL_RESOLUTION, state.data()));

            // Add contribution if the cell is fully inside
            if (result & 0b100) {
                return true;
            } else if (result & 0b10) {
                Value cellVolume = 1;
                for (unsigned i_dim=0; i_dim<Dimension; ++i_dim) {
                    cellVolume *= lengths[i_dim] * (CELL_RESOLUTION - 1);
                }
                volume += cellVolume;
            }
            return false;
        };

        tree.scan(visitor);
        if (!comparison.equal(volume, referenceVolume)) {
            throw std::runtime_error("Incorrect results " + std::to_string(volume) + " != " + std::to_string(referenceVolume));
        }
    }
}


template <class TGeometry>
void flexibleTree(benchmark::State& r_state)
{
    constexpr unsigned Dimension = TGeometry::Dimension;
    using Value = typename TGeometry::Coordinate;
    using Cell = cie::geo::Cell<TGeometry>;
    using Node = cie::geo::SpaceTreeNode<Cell,bool>;

    using Sampler = cie::geo::CartesianGridSampler<TGeometry>;
    using SplitPolicy = cie::geo::MidPointSplitPolicy<typename Node::SamplePointIterator,typename Node::value_iterator>;

    const cie::utils::Comparison<Value> comparison(
        std::numeric_limits<Value>::min(),
        5e-1
    );
    const Value referenceVolume =
        Dimension == 1 ? 2 :
        Dimension == 2 ? std::numbers::pi :
        Dimension == 3 ? 4.0 / 3.0 * std::numbers::pi :
        std::numeric_limits<Value>::max();

    const typename Node::Target target = [](const typename TGeometry::Point& r_sample) -> bool {
        return isInside(r_sample.data(), r_sample.data() + Dimension);
    };

    for ([[maybe_unused]] auto dummy : r_state) {
        // Construct tree root
        typename Sampler::Resolution resolution;
        std::fill(resolution.begin(), resolution.end(), CELL_RESOLUTION);
        Node root(typename Node::sampler_ptr(new Sampler(resolution)),
                  typename Node::split_policy_ptr(new SplitPolicy),
                  0,
                  makeRootGeometry<TGeometry>());

        // Build tree and compute volume
        std::atomic<Value> volume = 0;
        const std::function<void(const Node&)> visitor = [&volume, &target](const Node& r_node) -> void {
            if (r_node.isLeaf() && !r_node.isBoundary() && target(r_node.base())) {
                Value cellVolume = 1;
                if constexpr (cie::concepts::Cube<TGeometry>) {
                    cellVolume = std::pow(r_node.length(), Dimension);
                } else {
                    cellVolume = std::accumulate(r_node.lengths().begin(),
                                                 r_node.lengths().end(),
                                                 1,
                                                 std::multiplies<Value>());
                }
                volume += cellVolume;
            }
        };
        root.scan(target, visitor, MAX_DEPTH);
        if (!comparison.equal(volume, referenceVolume)) {
            throw std::runtime_error("Incorrect results " + std::to_string(volume) + " != " + std::to_string(referenceVolume));
        }
    }
}


BENCHMARK_TEMPLATE(contiguousTree, cie::geo::Cube<1,float>)->Name("ContiguousTree_float_1D_cube");
BENCHMARK_TEMPLATE(contiguousTree, cie::geo::Cube<2,float>)->Name("ContiguousTree_float_2D_cube");
BENCHMARK_TEMPLATE(contiguousTree, cie::geo::Cube<3,float>)->Name("ContiguousTree_float_3D_cube");

BENCHMARK_TEMPLATE(contiguousTree, cie::geo::Cube<1,double>)->Name("ContiguousTree_double_1D_cube");
BENCHMARK_TEMPLATE(contiguousTree, cie::geo::Cube<2,double>)->Name("ContiguousTree_double_2D_cube");
BENCHMARK_TEMPLATE(contiguousTree, cie::geo::Cube<3,double>)->Name("ContiguousTree_double_3D_cube");

BENCHMARK_TEMPLATE(contiguousTree, cie::geo::Box<1,float>)->Name("ContiguousTree_float_1D_box");
BENCHMARK_TEMPLATE(contiguousTree, cie::geo::Box<2,float>)->Name("ContiguousTree_float_2D_box");
BENCHMARK_TEMPLATE(contiguousTree, cie::geo::Box<3,float>)->Name("ContiguousTree_float_3D_box");

BENCHMARK_TEMPLATE(contiguousTree, cie::geo::Box<1,double>)->Name("ContiguousTree_double_1D_box");
BENCHMARK_TEMPLATE(contiguousTree, cie::geo::Box<2,double>)->Name("ContiguousTree_double_2D_box");
BENCHMARK_TEMPLATE(contiguousTree, cie::geo::Box<3,double>)->Name("ContiguousTree_double_3D_box");


BENCHMARK_TEMPLATE(flexibleTree, cie::geo::Cube<1,float>)->Name("Tree_float_1D_cube");
BENCHMARK_TEMPLATE(flexibleTree, cie::geo::Cube<2,float>)->Name("Tree_float_2D_cube");
BENCHMARK_TEMPLATE(flexibleTree, cie::geo::Cube<3,float>)->Name("Tree_float_3D_cube");

BENCHMARK_TEMPLATE(flexibleTree, cie::geo::Cube<1,double>)->Name("Tree_double_1D_cube");
BENCHMARK_TEMPLATE(flexibleTree, cie::geo::Cube<2,double>)->Name("Tree_double_2D_cube");
BENCHMARK_TEMPLATE(flexibleTree, cie::geo::Cube<3,double>)->Name("Tree_double_3D_cube");

BENCHMARK_TEMPLATE(flexibleTree, cie::geo::Box<1,float>)->Name("Tree_float_1D_box");
BENCHMARK_TEMPLATE(flexibleTree, cie::geo::Box<2,float>)->Name("Tree_float_2D_box");
BENCHMARK_TEMPLATE(flexibleTree, cie::geo::Box<3,float>)->Name("Tree_float_3D_box");

BENCHMARK_TEMPLATE(flexibleTree, cie::geo::Box<1,double>)->Name("Tree_double_1D_box");
BENCHMARK_TEMPLATE(flexibleTree, cie::geo::Box<2,double>)->Name("Tree_double_2D_box");
BENCHMARK_TEMPLATE(flexibleTree, cie::geo::Box<3,double>)->Name("Tree_double_3D_box");
