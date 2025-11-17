#ifndef CIE_FEM_ANSATZ_SPACE_HPP
#define CIE_FEM_ANSATZ_SPACE_HPP

// --- FEM Includes ---
#include "packages/maths/inc/Expression.hpp"
#include "packages/io/inc/GraphML.hpp"

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

    using typename ExpressionTraits<Value>::ConstSpan;

    using typename ExpressionTraits<Value>::Span;

public:
    AnsatzSpaceDerivative() noexcept = default;

    void evaluate(ConstSpan in, Span out) const;

    unsigned size() const noexcept;

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
class BufferedAnsatzSpace : public ExpressionTraits<typename TScalarExpression::Value>
{
private:
    using Base = ExpressionTraits<typename TScalarExpression::Value>;

public:
    static constexpr unsigned Dimension = Dim;

    using typename Base::Value;

    using typename Base::Span;

    using typename Base::ConstSpan;

    BufferedAnsatzSpace() noexcept;

    BufferedAnsatzSpace(std::span<const TScalarExpression> ansatzSet) noexcept;

    BufferedAnsatzSpace(std::span<const TScalarExpression> ansatzSet,
                        Span buffer);

    void evaluate(ConstSpan in, Span out) const;

    unsigned size() const noexcept;

    unsigned getMinBufferSize() const noexcept;

    void setBuffer(Span buffer);

    std::span<const TScalarExpression> ansatzSet() const noexcept;

private:
    std::span<unsigned,Dim> getIndexBuffer() const noexcept;

    Span getValueBuffer() const noexcept;

    std::span<const TScalarExpression> _set;

    mutable Span _buffer;
}; // class AnsatzSpace



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

    using typename Base::Span;

    using typename Base::ConstSpan;

    using AnsatzSet = DynamicArray<TScalarExpression>;

    using Derivative = AnsatzSpaceDerivative<TScalarExpression,Dim>;

public:
    AnsatzSpace() noexcept;

    AnsatzSpace(AnsatzSet&& rSet) noexcept;

    AnsatzSpace(const AnsatzSet& rSet);

    void evaluate(ConstSpan in, Span out) const;

    Derivative makeDerivative() const;

    unsigned size() const noexcept;

    std::span<const TScalarExpression> ansatzSet() const noexcept;

private:
    friend class AnsatzSpaceDerivative<TScalarExpression,Dim>;

    AnsatzSet _set;

    DynamicArray<Value> _buffer;

    BufferedAnsatzSpace<TScalarExpression,Dim> _wrapped;
}; // class AnsatzSpace


} // namespace cie::fem::maths



// --- IO --- //



namespace cie::fem::io {


template <class TScalarExpression, unsigned Dimension>
struct io::GraphML::Serializer<maths::AnsatzSpace<TScalarExpression,Dimension>>
{
    using Value = maths::AnsatzSpace<TScalarExpression,Dimension>;

    void header(Ref<XMLElement> rElement);

    void operator()(Ref<XMLElement> rElement,
                    Ref<const maths::AnsatzSpace<TScalarExpression,Dimension>> rInstance);
}; // struct GraphML::Serializer<AnsatzSpace>


template <class TScalarExpression, unsigned Dimension>
struct io::GraphML::Deserializer<maths::AnsatzSpace<TScalarExpression,Dimension>>
    : public io::GraphML::DeserializerBase<maths::AnsatzSpace<TScalarExpression,Dimension>>
{
    using Value = maths::AnsatzSpace<TScalarExpression,Dimension>;

    using io::GraphML::DeserializerBase<Value>::DeserializerBase;

    static void onElementBegin(Ptr<void> pThis,
                               std::string_view elementName,
                               std::span<GraphML::AttributePair> attributes);

    static void onText(Ptr<void> pThis,
                       std::string_view data);

    static void onElementEnd(Ptr<void> pThis,
                             std::string_view elementName);

private:
    typename Value::AnsatzSet _set;
}; // struct GraphML::Deserializer<AnsatzSpace>


} // namespace cie::fem::io


#endif

#include "packages/maths/impl/AnsatzSpace_impl.hpp"
