#pragma once

// --- Utility Includes ---
#include "packages/stl_extension/inc/DynamicArray.hpp"
#include "packages/types/inc/types.hpp"

// --- STL Includes ---
#include <utility>


namespace cie::fem {


///@addtogroup fem
///@{

template <concepts::Numeric NT>
class QuadratureBase {
public:
    using Value = NT;

    using NodeContainer   = DynamicArray<NT>;

    using WeightContainer = DynamicArray<NT>;

public:
    Size numberOfNodes() const;

    Size order() const;

    Ref<const NodeContainer> nodes() const;

    Ref<const WeightContainer> weights() const;

protected:
    QuadratureBase() noexcept = default;

    QuadratureBase(NodeContainer&& rNodes,
                   WeightContainer&& rWeights);

    QuadratureBase(std::pair<NodeContainer, WeightContainer>&& rNodesAndWeights);

    QuadratureBase(Ref<const NodeContainer> rNodes,
                   Ref<const WeightContainer> rWeights);

private:
    NodeContainer   _nodes;

    WeightContainer _weights;
};

///@}

} // namespace cie::fem
