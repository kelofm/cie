// --- Utility Includes ---
#include "packages/macros/inc/exceptions.hpp"

// --- FEM Includes ---
#include "packages/utilities/inc/template_macros.hpp"

// --- Internal Includes ---
#include "packages/numeric/inc/QuadratureBase.hpp"


namespace cie::fem {


template <concepts::Numeric NT>
QuadratureBase<NT>::QuadratureBase(typename QuadratureBase<NT>::NodeContainer&& rNodes,
                                   typename QuadratureBase<NT>::WeightContainer&& rWeights)
    : _nodes( std::move(rNodes) ),
      _weights( std::move(rWeights) )
{
}


template <concepts::Numeric NT>
QuadratureBase<NT>::QuadratureBase(std::pair<typename QuadratureBase<NT>::NodeContainer,
                                             typename QuadratureBase<NT>::WeightContainer>&& rNodesAndWeights )
    : _nodes( std::move(rNodesAndWeights.first) ),
      _weights( std::move(rNodesAndWeights.second) )
{
}


template <concepts::Numeric NT>
QuadratureBase<NT>::QuadratureBase(Ref<const typename QuadratureBase<NT>::NodeContainer> rNodes,
                                   Ref<const typename QuadratureBase<NT>::WeightContainer> rWeights)
    : _nodes( rNodes ),
      _weights( rWeights )
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
