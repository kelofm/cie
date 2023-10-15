#ifndef CIE_UTILS_MATHS_POLYNOMIAL_HPP
#define CIE_UTILS_MATHS_POLYNOMIAL_HPP

// --- Utility Includes ---
#include "packages/compile_time/packages/concepts/inc/iterator_concepts.hpp"


namespace cie::utils {


/** @brief Evaluate a polynomial as it is written in its canonical form.
 *  @tparam TValue: value type of the argument and the result.
 *  @tparam TItBegin: coefficient iterator type.
 *  @tparam TItEnd: coefficient iterator type.
 *  @param argument: position at which to evaluate the polynomial.
 *  @param it_begin: iterator pointing to the coefficient of the lowest power (0).
 *  @param it_end: iterator past the last coefficient.
 *  @note The coefficients are expected to be sorted with respect to their order in the polynomial
 *        (increasing order, starting with order 0 and including null value coefficients).
 *  @warning This function does not check for overflows.
 *  @warning This scheme is prone to excessive numerical errors for higher order terms
 *           or large arguments.
 *  @ingroup cieutils
 */
template <concepts::Numeric TValue, concepts::WeakIterator<TValue> TItBegin, concepts::WeakIterator<TValue> TItEnd>
TValue evaluatePolynomialNaive(TValue argument, TItBegin it_begin, const TItEnd it_end) noexcept;


/** @brief Evaluate a polynomial as it is written in its canonical form.
 *  @tparam TValue: value type of the argument and the result.
 *  @tparam TItBegin: coefficient iterator type.
 *  @tparam TItEnd: coefficient iterator type.
 *  @param argument: position at which to evaluate the polynomial.
 *  @param it_begin: iterator pointing to the coefficient of the lowest power (0).
 *  @param it_end: iterator past the last coefficient.
 *  @note The coefficients are expected to be sorted with respect to their order in the polynomial
 *        (increasing order, starting with order 0 and including null value coefficients).
 *  @warning This function does not check for overflows.
 *  @warning This scheme can become inaccurate outside the [-1,1] argument range.
 *           Use @ref evaluatePolynomialHornerStabilized in these cases instead.
 *  @ingroup cieutils
 */
template <concepts::Numeric TValue, concepts::WeakIterator<TValue> TItBegin, concepts::WeakIterator<TValue> TItEnd>
TValue evaluatePolynomialHorner(TValue argument, const TItBegin it_begin, TItEnd it_end) noexcept;


/** @brief Evaluate a polynomial as it is written in its canonical form.
 *  @tparam TValue: value type of the argument and the result.
 *  @tparam TItBegin: coefficient iterator type.
 *  @tparam TItEnd: coefficient iterator type.
 *  @param argument: position at which to evaluate the polynomial.
 *  @param it_begin: iterator pointing to the coefficient of the lowest power (0).
 *  @param it_end: iterator past the last coefficient.
 *  @note The coefficients are expected to be sorted with respect to their order in the polynomial
 *        (increasing order, starting with order 0 and including null value coefficients).
 *  @note This scheme sacrifices some performance for accuracy outside the [-1,1] argument range.
 *        Use @ref evaluatePolynomialHorner if your arguments are normalized to the [-1,1] range.
 *  @warning This function does not check for overflows.
 *  @warning This function computes the inverse of @a argument but cuts off the operation
 *           is below machine epsilon, returning 0.
 *  @ingroup cieutils
 */
template <concepts::Numeric TValue, concepts::WeakIterator<TValue> TItBegin, concepts::WeakIterator<TValue> TItEnd>
TValue evaluatePolynomialHornerStabilized(TValue argument, TItBegin it_begin, const TItEnd it_end) noexcept;


/** @brief Evaluate a polynomial as it is written in its canonical form.
 *  @tparam TValue: value type of the argument and the result.
 *  @tparam TItBegin: coefficient iterator type.
 *  @tparam TItEnd: coefficient iterator type.
 *  @param argument: position at which to evaluate the polynomial.
 *  @param it_begin: iterator pointing to the coefficient of the lowest power (0).
 *  @param it_end: iterator past the last coefficient.
 *  @note The coefficients are expected to be sorted with respect to their order in the polynomial
 *        (increasing order, starting with order 0 and including null value coefficients).
 *  @note This scheme has a lower performance outside the [-1,1] argument range.
 *        Use @ref evaluatePolynomialHorner if your arguments are normalized to the [-1,1] range.
 *  @warning This function does not check for overflows.
 *  @ingroup cieutils
 */
template <concepts::Numeric TValue, concepts::WeakIterator<TValue> TItBegin, concepts::WeakIterator<TValue> TItEnd>
TValue evaluatePolynomialHornerCompound(TValue argument, TItBegin it_begin, TItEnd it_end) noexcept;


} // namespace cie::utils

#include "packages/maths/impl/polynomial_evaluation_impl.hpp"


#endif
