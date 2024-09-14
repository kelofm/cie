#ifndef CIE_FEM_ANSATZ_SPACE_HPP
#define CIE_FEM_ANSATZ_SPACE_HPP

// --- FEM Includes ---
#include "packages/maths/inc/Expression.hpp"

// --- Utility Includes ---
#include "packages/stl_extension/inc/DynamicArray.hpp"
#include "packages/stl_extension/inc/StaticArray.hpp"
#include "packages/concurrency/inc/ThreadLocal.hpp"


namespace cie::fem::maths {


template <class TScalarExpression, unsigned Dim>
class AnsatzSpace;



template <class TScalarExpression, unsigned Dim>
class AnsatzSpaceDerivative : public ExpressionTraits<typename TScalarExpression::Value>
{
public:
    static constexpr unsigned Dimension = Dim;

    using typename ExpressionTraits<typename TScalarExpression::Value>::Value;

    using typename ExpressionTraits<Value>::ConstIterator;

    using typename ExpressionTraits<Value>::Iterator;

public:
    AnsatzSpaceDerivative() noexcept = default;

    void evaluate(ConstIterator itArgumentBegin,
                  ConstIterator itArgumentEnd,
                  Iterator itOut) const;

    Size size() const noexcept;

private:
    friend class AnsatzSpace<TScalarExpression,Dim>;

    AnsatzSpaceDerivative(Ref<const AnsatzSpace<TScalarExpression,Dim>> rAnsatzSpace);

private:
    DynamicArray<TScalarExpression> _ansatzSet;

    DynamicArray<typename TScalarExpression::Derivative> _derivativeSet;

    using IndexBuffer = StaticArray<unsigned,Dim>;

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
template <class TScalarExpression, unsigned Dim>
class AnsatzSpace : public ExpressionTraits<typename TScalarExpression::Value>
{
private:
    using Base = ExpressionTraits<typename TScalarExpression::Value>;

public:
    static constexpr unsigned Dimension = Dim;

    using typename Base::Value;

    using typename Base::Iterator;

    using typename Base::ConstIterator;

    using AnsatzSet = DynamicArray<TScalarExpression>;

    using Derivative = AnsatzSpaceDerivative<TScalarExpression,Dim>;

public:
    AnsatzSpace() noexcept;

    AnsatzSpace(AnsatzSet&& rSet) noexcept;

    AnsatzSpace(const AnsatzSet& rSet);

    void evaluate(ConstIterator itArgumentBegin,
                  ConstIterator itArgumentEnd,
                  Iterator itOut) const;

    Derivative makeDerivative() const;

    unsigned size() const noexcept;

private:
    AnsatzSet _set;

    using IndexBuffer = StaticArray<unsigned,Dim>;

    using ValueBuffer = DynamicArray<Value>;

    friend class AnsatzSpaceDerivative<TScalarExpression,Dim>;

    /// @brief A threadsafe container for eliminating allocations from @ref AnsatzSpace::evaluate.
    mutable mp::ThreadLocal<
        IndexBuffer,
        ValueBuffer
    > _buffer;
}; // class AnsatzSpace


} // namespace cie::fem::maths


#endif

#include "packages/maths/impl/AnsatzSpace_impl.hpp"
