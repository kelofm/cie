// --- Utility Includes ---
#include "packages/macros/inc/exceptions.hpp"

// --- FEM Includes ---
#include "packages/utilities/inc/template_macros.hpp"

// --- Internal Includes ---
#include "packages/numeric/inc/QuadratureBase.hpp"


namespace cie::fem {


template <concepts::Numeric NT>
QuadratureBase<NT>::QuadratureBase(typename QuadratureBase<NT>::NodeContainer&& r_nodes,
                                   typename QuadratureBase<NT>::WeightContainer&& r_weights)
    : _nodes( std::move(r_nodes) ),
      _weights( std::move(r_weights) )
{
}


template <concepts::Numeric NT>
QuadratureBase<NT>::QuadratureBase(std::pair<typename QuadratureBase<NT>::NodeContainer,
                                             typename QuadratureBase<NT>::WeightContainer>&& r_nodesAndWeights )
    : _nodes( std::move(r_nodesAndWeights.first) ),
      _weights( std::move(r_nodesAndWeights.second) )
{
}


template <concepts::Numeric NT>
QuadratureBase<NT>::QuadratureBase(Ref<const typename QuadratureBase<NT>::NodeContainer> r_nodes,
                                   Ref<const typename QuadratureBase<NT>::WeightContainer> r_weights)
    : _nodes( r_nodes ),
      _weights( r_weights )
{
}


template <concepts::Numeric NT>
Size QuadratureBase<NT>::numberOfNodes() const
{
    return _nodes.size();
}


template <concepts::Numeric NT>
Size QuadratureBase<NT>::order() const
{
    return _nodes.size();
}


template <concepts::Numeric NT>
Ref<const typename QuadratureBase<NT>::NodeContainer>
QuadratureBase<NT>::nodes() const
{
    return _nodes;
}


template <concepts::Numeric NT>
Ref<const typename QuadratureBase<NT>::WeightContainer>
QuadratureBase<NT>::weights() const
{
    return _weights;
}


CIE_FEM_INSTANTIATE_NUMERIC_TEMPLATE(QuadratureBase);


} // namespace cie::fem
