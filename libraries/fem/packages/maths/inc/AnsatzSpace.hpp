#ifndef CIE_FEM_ANSATZ_SPACE_HPP
#define CIE_FEM_ANSATZ_SPACE_HPP

// --- FEM Includes ---
#include "packages/maths/inc/Expression.hpp"

// --- Utility Includes ---
#include "packages/stl_extension/inc/DynamicArray.hpp"
#include "packages/stl_extension/inc/StaticArray.hpp"
#include "packages/concurrency/inc/ThreadLocal.hpp"


namespace cie::fem::maths {


template <class TScalarExpression, unsigned Dimension>
class AnsatzSpace;



template <class TScalarExpression, unsigned Dimension>
class AnsatzSpaceDerivative : public ExpressionTraits<typename TScalarExpression::Value>
{
public:
    using typename ExpressionTraits<typename TScalarExpression::Value>::Value;

    using typename ExpressionTraits<Value>::ConstIterator;

    using typename ExpressionTraits<Value>::Iterator;

public:
    AnsatzSpaceDerivative() noexcept = default;

    void evaluate(ConstIterator it_argumentBegin,
                  ConstIterator it_argumentEnd,
                  Iterator it_out) const;

    Size size() const noexcept;

private:
    friend class AnsatzSpace<TScalarExpression,Dimension>;

    AnsatzSpaceDerivative(Ref<const AnsatzSpace<TScalarExpression,Dimension>> r_ansatzSpace);

private:
    DynamicArray<TScalarExpression> _ansatzSet;

    DynamicArray<typename TScalarExpression::Derivative> _derivativeSet;

    using IndexBuffer = StaticArray<unsigned,Dimension>;

    using ValueBuffer = DynamicArray<Value>;

    /// @brief A threadsafe container for eliminating allocations from @ref AnsatzSpaceDerivative::evaluate.
    mutable mp::ThreadLocal<
        IndexBuffer, // <== indices for the cartesian product
        ValueBuffer, // <== buffer for ansatz function values at the cartesian grid points
        ValueBuffer  // <== buffer for the derivatives at the cartesian grid points
    > _buffer;
}; // class AnsatzSpaceDerivative



/** @brief A set of multidimensional functions constructed from the cartesian product of a set of scalar basis functions.
 */
template <class TScalarExpression, unsigned Dimension>
class AnsatzSpace : public ExpressionTraits<typename TScalarExpression::Value>
{
private:
    using Base = ExpressionTraits<typename TScalarExpression::Value>;

public:
    using typename Base::Value;

    using typename Base::Iterator;

    using typename Base::ConstIterator;

    using AnsatzSet = DynamicArray<TScalarExpression>;

public:
    AnsatzSpace() noexcept;

    AnsatzSpace(AnsatzSet&& r_set) noexcept;

    AnsatzSpace(const AnsatzSet& r_set);

    void evaluate(ConstIterator it_argumentBegin,
                  ConstIterator it_argumentEnd,
                  Iterator it_out) const;

    AnsatzSpaceDerivative<TScalarExpression,Dimension> makeDerivative() const;

    unsigned size() const noexcept;

private:
    AnsatzSet _set;

    using IndexBuffer = StaticArray<unsigned,Dimension>;

    using ValueBuffer = DynamicArray<Value>;

    friend class AnsatzSpaceDerivative<TScalarExpression,Dimension>;

    /// @brief A threadsafe container for eliminating allocations from @ref AnsatzSpace::evaluate.
    mutable mp::ThreadLocal<
        IndexBuffer,
        ValueBuffer
    > _buffer;
}; // class AnsatzSpace


} // namespace cie::fem::maths


#endif

#include "packages/maths/impl/AnsatzSpace_impl.hpp"
