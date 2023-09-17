#ifndef CIE_FEM_ANSATZ_SPACE_IMPL_HPP
#define CIE_FEM_ANSATZ_SPACE_IMPL_HPP

// --- FEM Includes ---
#include "packages/maths/inc/AnsatzSpace.hpp"
#include "packages/maths/inc/CartesianProduct.hpp"

// --- Utility Includes ---
#include "packages/macros/inc/checks.hpp"
#include "packages/maths/inc/power.hpp"
#include "packages/stl_extension/inc/StaticArray.hpp"

// --- STL Includes ---
#include <algorithm>


namespace cie::fem::maths {


template <class TScalarExpression, unsigned Dimension>
inline void AnsatzSpaceDerivative<TScalarExpression,Dimension>::evaluate(ConstIterator it_argumentBegin,
                                                                         ConstIterator it_argumentEnd,
                                                                         Iterator it_out) const
{
    const unsigned setSize = _ansatzSet.size();
    CIE_OUT_OF_RANGE_CHECK(setSize == _derivativeSet.size())
    CIE_OUT_OF_RANGE_CHECK(std::distance(it_argumentBegin, it_argumentEnd) == Dimension)

    Ref<IndexBuffer> r_indexBuffer      = _buffer.template get<0>();
    Ref<ValueBuffer> r_valueBuffer      = _buffer.template get<1>();
    Ref<ValueBuffer> r_derivativeBuffer = _buffer.template get<2>();

    // Fill the value and derivative buffers
    {
        Ptr<Value> p_value      = r_valueBuffer.data();
        Ptr<Value> p_derivative = r_derivativeBuffer.data();

        for (auto it=it_argumentBegin; it!=it_argumentEnd; ++it) {
            const auto it_end = it + 1;
            for (const auto& r_scalarExpression : _ansatzSet) {
                r_scalarExpression.evaluate(it,
                                            it_end,
                                            p_value++);
            } // for scalarExpression in ansatzSet
            for (const auto& r_scalarExpression : _derivativeSet) {
                r_scalarExpression.evaluate(it,
                                            it_end,
                                            p_derivative++);
            } // for scalarExpression in derivativeSet
        } // for component in arguments
    } // fill the value and derivative buffers

    // Compute the modified cartesian product
    for (unsigned i_derivative=0; i_derivative<Dimension; ++i_derivative) {
        do {
            *it_out = static_cast<Value>(1);
            unsigned i_index = 0;

            // First loop through the value buffer until
            // the derivative index is hit
            for (; i_index<i_derivative; ++i_index) {
                *it_out *= r_valueBuffer[r_indexBuffer[i_index] + i_index * setSize];
            }

            // Then, use the derivative
            *it_out *= r_derivativeBuffer[r_indexBuffer[i_index] + i_index * setSize];

            // Finally, loop through the rest
            // of the value buffer
            for (++i_index; i_index<r_indexBuffer.size(); ++i_index) {
                *it_out *= r_valueBuffer[r_indexBuffer[i_index] + i_index * setSize];
            }

            ++it_out;
        } while (CartesianProduct<Dimension>::next(setSize, r_indexBuffer.begin()));
    } // for i_derivative in range(Dimension)
}


template <class TScalarExpression, unsigned Dimension>
Size AnsatzSpaceDerivative<TScalarExpression,Dimension>::size() const noexcept
{
    return intPow(_ansatzSet.size(), Dimension) * Dimension;
}


template <class TScalarExpression, unsigned Dimension>
AnsatzSpaceDerivative<TScalarExpression,Dimension>::AnsatzSpaceDerivative(Ref<const AnsatzSpace<TScalarExpression,Dimension>> r_ansatzSpace)
    : _ansatzSet(r_ansatzSpace._set),
      _derivativeSet(),
      _buffer(IndexBuffer(),
              ValueBuffer(intPow(_ansatzSet.size(), Dimension)),
              ValueBuffer(intPow(_ansatzSet.size(), Dimension)))
{
    Ref<IndexBuffer> r_indexBuffer = _buffer.template get<0>();
    std::fill(r_indexBuffer.begin(),
              r_indexBuffer.end(),
              0u);

    _derivativeSet.reserve(_ansatzSet.size());
    for (const auto& r_scalarExpression : _ansatzSet) {
        _derivativeSet.emplace_back(r_scalarExpression.makeDerivative());
    }
}


template <class TScalarExpression, unsigned Dimension>
AnsatzSpace<TScalarExpression,Dimension>::AnsatzSpace() noexcept
    : AnsatzSpace({})
{
}


template <class TScalarExpression, unsigned Dimension>
AnsatzSpace<TScalarExpression,Dimension>::AnsatzSpace(AnsatzSet&& r_set) noexcept
    : _set(std::move(r_set)),
      _buffer(IndexBuffer(),
              ValueBuffer(intPow(_set.size(), Dimension)))
{
    std::fill(_buffer.template get<0>().begin(),
              _buffer.template get<0>().end(),
              0u);
}


template <class TScalarExpression, unsigned Dimension>
AnsatzSpace<TScalarExpression,Dimension>::AnsatzSpace(const AnsatzSet& r_set)
    : AnsatzSpace(AnsatzSet(r_set))
{
}


template <class TScalarExpression, unsigned Dimension>
inline void AnsatzSpace<TScalarExpression,Dimension>::evaluate(ConstIterator it_argumentBegin,
                                                               ConstIterator it_argumentEnd,
                                                               Iterator it_out) const
{
    CIE_OUT_OF_RANGE_CHECK(std::distance(it_argumentBegin, it_argumentEnd) == Dimension)
    Ref<IndexBuffer> r_indexBuffer = _buffer.template get<0>();
    Ref<ValueBuffer> r_valueBuffer = _buffer.template get<1>();

    // No need to clear the index buffer of its leftover state
    // (the leftover state is all zeros)
    //std::fill(r_indexBuffer->begin(), r_indexBuffer->end(), 0);

    // Fill the value buffer
    Ptr<Value> p_value = r_valueBuffer.data();
    for (; it_argumentBegin!=it_argumentEnd; it_argumentBegin++) {
        for (const auto& r_scalarExpression : _set) {
            r_scalarExpression.evaluate(it_argumentBegin,
                                        it_argumentBegin + 1,
                                        p_value++);
        } // for expression in expressions
    } // for argument in arguments

    const unsigned setSize = _set.size();
    do {
        // Compute product of bases
        *it_out = static_cast<Value>(1);
        for (unsigned i_index=0; i_index<r_indexBuffer.size(); ++i_index) {
            *it_out *= r_valueBuffer.at(r_indexBuffer.at(i_index) + i_index * setSize);
        }
        ++it_out;
    } while (CartesianProduct<Dimension>::next(setSize, r_indexBuffer.data()));
}


template <class TScalarExpression, unsigned Dimension>
typename AnsatzSpace<TScalarExpression,Dimension>::Derivative
AnsatzSpace<TScalarExpression,Dimension>::makeDerivative() const
{
    return Derivative(*this);
}


template <class TScalarExpression, unsigned Dimension>
unsigned AnsatzSpace<TScalarExpression,Dimension>::size() const noexcept
{
    return intPow(_set.size(), Dimension);
}


} // namespace cie::fem::maths



#endif
