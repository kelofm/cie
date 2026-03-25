#pragma once

// --- FEM Includes ---
#include "packages/maths/inc/Expression.hpp"


namespace cie::fem {


template <maths::Expression TIntegrand, maths::JacobianExpression TJacobian>
class TransformedIntegrand : public maths::ExpressionTraits<typename TIntegrand::Value> {
public:
    static constexpr bool IsStatic = maths::StaticExpression<TIntegrand> && maths::StaticExpression<TJacobian>;

    using typename maths::ExpressionTraits<typename TIntegrand::Value>::Value;

    using typename maths::ExpressionTraits<typename TIntegrand::Value>::ConstSpan;

    using typename maths::ExpressionTraits<typename TIntegrand::Value>::Span;

    using typename maths::ExpressionTraits<typename TIntegrand::Value>::BufferSpan;

    using Integrand = TIntegrand;

    using Jacobian = TJacobian;

    TransformedIntegrand() noexcept {}

    TransformedIntegrand(
        RightRef<TIntegrand> rIntegrand,
        RightRef<TJacobian> rInverseJacobian) noexcept
            : _integrand(std::move(rIntegrand)),
              _inverseJacobian(std::move(rInverseJacobian))
    {}

    unsigned size() const noexcept
    requires (!IsStatic) {
        return _integrand.size();}

    static constexpr unsigned size() noexcept
    requires (IsStatic) {
        return TIntegrand::size();}

    unsigned bufferSize() const noexcept
    requires (!IsStatic) {
        return _integrand.bufferSize() + _inverseJacobian.bufferSize();}

    static constexpr unsigned bufferSize() noexcept
    requires (IsStatic) {
        return TIntegrand::bufferSize() + Jacobian::bufferSize();}

    void evaluate(
        ConstSpan in,
        Span out,
        [[maybe_unused]] BufferSpan buffer) const {
            BufferSpan integrandBuffer, jacobianBuffer;
            [[maybe_unused]] std::array<Value,maths::StaticExpressionSize<TIntegrand>::bufferSize> integrandBufferArray;
            [[maybe_unused]] std::array<Value,maths::StaticExpressionSize<TJacobian>::bufferSize> jacobianBufferArray;

            if constexpr (IsStatic) {
                integrandBuffer = integrandBufferArray;
                jacobianBuffer  = jacobianBufferArray;
            } else {
                const unsigned integrandBufferSize = _integrand.bufferSize();
                integrandBuffer = BufferSpan(
                    buffer.data(),
                    integrandBufferSize);
                jacobianBuffer = BufferSpan(
                    buffer.data() + integrandBufferSize,
                    buffer.data() + buffer.size());
            }

            const Value determinant = std::abs(_inverseJacobian.evaluateDeterminant(
                in,
                jacobianBuffer));
            //const Value scale = static_cast<Value>(1) / determinant;
            const Value scale = determinant;
            _integrand.evaluate(
                in,
                out,
                integrandBuffer);
            for (Value& rComponent : out) {
                rComponent *= scale;
            }
    }

private:
    TIntegrand _integrand;

    TJacobian _inverseJacobian;
}; // class TransformedIntegrand


template <maths::Expression TIntegrand, maths::JacobianExpression TJacobian>
TransformedIntegrand<TIntegrand,TJacobian>
makeTransformedIntegrand(
    RightRef<TIntegrand> rIntegrand,
    RightRef<TJacobian> rJacobian) {
        return TransformedIntegrand(
            std::move(rIntegrand),
            std::move(rJacobian));
}


} // namespace cie::fem
