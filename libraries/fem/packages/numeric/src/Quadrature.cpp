// --- FEM Includes ---
#include "packages/numeric/inc/Quadrature.hpp"
#include "packages/maths/inc/OuterProduct.hpp"
#include "packages/utilities/inc/template_macros.hpp"

// --- Utility Includes ---
#include "packages/maths/inc/power.hpp"
#include "packages/macros/inc/checks.hpp"

// --- STL Includes ---
#include <algorithm>


namespace cie::fem {


template <concepts::Numeric TValue, unsigned Dimension>
Quadrature<TValue,Dimension>::Quadrature(Ref<const QuadratureBase<TValue>> r_base)
    : Quadrature(r_base.nodes(), r_base.weights())
{
}


template <concepts::Numeric TValue, unsigned Dimension>
Quadrature<TValue,Dimension>::Quadrature(Ref<const typename QuadratureBase<TValue>::NodeContainer> r_nodes,
                                         Ref<const typename QuadratureBase<TValue>::WeightContainer> r_weights)
    : _nodesAndWeights()
{
    const unsigned numberOfNodes = r_nodes.size();
    CIE_OUT_OF_RANGE_CHECK(numberOfNodes == r_weights.size())

    if (!r_nodes.empty()) {
        // Resize nodes and weights
        this->_nodesAndWeights.resize(intPow(r_nodes.size(), Dimension));

        // Create an index buffer for constructing the outer product
        StaticArray<unsigned,Dimension> indexBuffer;
        std::fill(indexBuffer.begin(),
                  indexBuffer.end(),
                  0);

        // Construct the outer product
        auto it_item = this->_nodesAndWeights.begin();
        do {
            // Compute the current component of the outer product
            std::fill(it_item->begin(),
                      it_item->end(),
                      static_cast<TValue>(1));
            Ref<TValue> r_weight = it_item->back();
            for (unsigned i_index=0; i_index<indexBuffer.size(); ++i_index) {
                const auto i = indexBuffer[i_index];
                it_item->at(i_index) = r_nodes[i];
                r_weight *= r_weights[i];
            } // for index in indexBuffer
            ++it_item;
        } while (maths::OuterProduct<Dimension>::next(numberOfNodes, indexBuffer.data()));
    } // if r_nodes
}


CIE_FEM_INSTANTIATE_TEMPLATE(Quadrature);


} // namespace cie::fem
