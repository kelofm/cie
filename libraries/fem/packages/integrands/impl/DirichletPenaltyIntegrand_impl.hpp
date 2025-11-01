#pragma once

// --- External Includes ---
#include <Eigen/Dense> // Eigen::Map

// help the language server
#include "packages/integrands/inc/DirichletPenaltyIntegrand.hpp"


namespace cie::fem {


template <maths::Expression TDirichlet, maths::Expression TAnsatz>
DirichletPenaltyIntegrand<TDirichlet,TAnsatz>::DirichletPenaltyIntegrand()
    : DirichletPenaltyIntegrand(TDirichlet(), 0, nullptr)
{
}


template <maths::Expression TDirichlet, maths::Expression TAnsatz>
DirichletPenaltyIntegrand<TDirichlet,TAnsatz>::DirichletPenaltyIntegrand(Ref<const TDirichlet> rDirichletFunctor,
                                                                         const Value penalty,
                                                                         Ref<const TAnsatz> rAnsatzSpace)
    : _penalty(penalty),
      _pDirichletFunctor(&rDirichletFunctor),
      _pAnsatzSpace(&rAnsatzSpace),
      _buffer()
{
}


template <maths::Expression TDirichlet, maths::Expression TAnsatz>
DirichletPenaltyIntegrand<TDirichlet,TAnsatz>::DirichletPenaltyIntegrand(Ref<const TDirichlet> rDirichletFunctor,
                                                                         const Value penalty,
                                                                         Ref<const TAnsatz> rAnsatzSpace,
                                                                         std::span<Value> buffer)
    : DirichletPenaltyIntegrand(rDirichletFunctor, penalty, rAnsatzSpace)
{
    this->setBuffer(buffer);
}


template <maths::Expression TDirichlet, maths::Expression TAnsatz>
void DirichletPenaltyIntegrand<TDirichlet,TAnsatz>::evaluate(ConstSpan in, Span out) const
{
    CIE_OUT_OF_RANGE_CHECK(this->getMinBufferSize() <= _buffer.size())
    CIE_CHECK_POINTER(_pAnsatzSpace)

    Ref<const TAnsatz> rAnsatzSpace = *_pAnsatzSpace;
    const unsigned ansatzCount = rAnsatzSpace.size();

    using EigenDenseMatrix = Eigen::Matrix<Value,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>;
    using EigenAdaptor = Eigen::Map<EigenDenseMatrix>;

    // Compute LHS contribution.
    Span ansatzBuffer(_buffer.data(), ansatzCount);
    rAnsatzSpace.evaluate(in, ansatzBuffer);
    EigenAdaptor ansatzAdaptor(_buffer.data(), ansatzCount, 1);
    EigenAdaptor lhsAdaptor(out.data(), ansatzCount, ansatzCount);

    lhsAdaptor = ansatzAdaptor * _penalty * ansatzAdaptor.transpose();

    // Compute RHS contribution.
    Span dirichletBuffer(_buffer.data() + ansatzCount, _pDirichletFunctor->size());
    _pDirichletFunctor->evaluate(in, dirichletBuffer);

    for (unsigned iStateComponent=0u; iStateComponent<dirichletBuffer.size(); ++iStateComponent) {
        const Value dirichlet = dirichletBuffer[iStateComponent];
        std::transform(ansatzBuffer.begin(),
                       ansatzBuffer.end(),
                       out.begin() + ansatzCount * ansatzCount + iStateComponent * ansatzCount,
                       [this, dirichlet] (const Value ansatz) -> Value {
                            return dirichlet * ansatz * _penalty;
                       });
    }
}


template <maths::Expression TDirichlet, maths::Expression TAnsatz>
unsigned DirichletPenaltyIntegrand<TDirichlet,TAnsatz>::size() const
{
    const auto ansatzCount = _pAnsatzSpace->size();
    return
          ansatzCount * ansatzCount     // <= LHS contribution
        + ansatzCount;                  // <= RHS contribution
}


template <maths::Expression TDirichlet, maths::Expression TAnsatz>
void DirichletPenaltyIntegrand<TDirichlet,TAnsatz>::setBuffer(std::span<Value> buffer)
{
    CIE_OUT_OF_RANGE_CHECK(this->getMinBufferSize() <= buffer.size())
    _buffer = buffer;
}


template <maths::Expression TDirichlet, maths::Expression TAnsatz>
unsigned DirichletPenaltyIntegrand<TDirichlet,TAnsatz>::getMinBufferSize() const noexcept
{
    return _pAnsatzSpace->size() + _pDirichletFunctor->size();
}


} // namespace cie::fem
