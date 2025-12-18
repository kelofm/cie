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
Quadrature<TValue,Dimension>::Quadrature(Ref<const QuadratureBase<TValue>> rBase)
    : Quadrature(rBase.nodes(), rBase.weights())
{
}


template <concepts::Numeric TValue, unsigned Dimension>
Quadrature<TValue,Dimension>::Quadrature(Ref<const typename QuadratureBase<TValue>::NodeContainer> rNodes,
                                         Ref<const typename QuadratureBase<TValue>::WeightContainer> rWeights)
    : _nodesAndWeights()
{
    const unsigned numberOfNodes = rNodes.size();
    CIE_OUT_OF_RANGE_CHECK(numberOfNodes == rWeights.size())

    if (!rNodes.empty()) {
        // Resize nodes and weights
        this->_nodesAndWeights.resize(intPow(rNodes.size(), Dimension));

        // Create an index buffer for constructing the outer product
        StaticArray<unsigned,Dimension> indexBuffer;
        std::fill_n(
            indexBuffer.data(),
            Dimension,
            0u);

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
                it_item->at(i_index) = rNodes[i];
                r_weight *= rWeights[i];
            } // for index in indexBuffer
            ++it_item;
        } while (cie::maths::OuterProduct<Dimension>::next(numberOfNodes, indexBuffer.data()));
    } // if rNodes
}


CIE_FEM_INSTANTIATE_TEMPLATE(Quadrature);


} // namespace cie::fem
