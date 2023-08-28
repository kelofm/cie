#ifndef CIE_FEM_NUMERIC_QUADRATURE_BASE_HPP
#define CIE_FEM_NUMERIC_QUADRATURE_BASE_HPP

// --- Utility Includes ---
#include "packages/stl_extension/inc/DynamicArray.hpp"
#include "packages/types/inc/types.hpp"

// --- STL Includes ---
#include <utility>


namespace cie::fem {


///@addtogroup fem
///@{

template <concepts::Numeric NT>
class QuadratureBase
{
public:
    using NodeContainer   = DynamicArray<NT>;

    using WeightContainer = DynamicArray<NT>;

public:
    Size numberOfNodes() const;

    Size order() const;

    Ref<const NodeContainer> nodes() const;

    Ref<const WeightContainer> weights() const;

protected:
    QuadratureBase() noexcept = default;

    QuadratureBase(NodeContainer&& r_nodes,
                   WeightContainer&& r_weights);

    QuadratureBase(std::pair<NodeContainer, WeightContainer>&& r_nodesAndWeights);

    QuadratureBase(Ref<const NodeContainer> r_nodes,
                   Ref<const WeightContainer> r_weights);

private:
    NodeContainer   _nodes;

    WeightContainer _weights;
};

///@}

} // namespace cie::fem

#endif