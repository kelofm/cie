#ifndef CIE_FEM_GRAPH_CONNECTIVITY_HPP
#define CIE_FEM_GRAPH_CONNECTIVITY_HPP

// --- FEM Includes ---
#include "packages/maths/inc/AnsatzSpace.hpp"
#include "packages/graph/inc/BoundaryID.hpp"
#include "packages/utilities/inc/kernel.hpp"

// --- Utility Includes ---
#include "packages/compile_time/packages/concepts/inc/functional.hpp"
#include "packages/stl_extension/inc/DynamicArray.hpp"
#include "packages/maths/inc/Comparison.hpp"

// --- STL Includes ---
#include <span>


namespace cie::fem {


/** @brief Collect ansatz functions that don't vanish on boundaries.
 *  @details Each boundary is scanned at sample points constructed from the cartesian
 *           product of the input sample points, and all ansatz functions that have at
 *           least one non-zero value at one of those points are considered to require
 *           connectivity on that boundary. An input functor is called with each
 *           @ref BoundaryID - ansatz function index pair exactly once.
 *  @param r_ansatzSpace @ref AnsatzSpace to scan the functions of.
 *  @param r_functor Functor that gets called with each @ref BoundaryID and non-vanishing
 *                   ansatz function index.
 *  @param p_sampleBegin Ptr to the beginning of the array of sample nodes to evaluate
 *                       the ansatz functions at.
 *  @param p_sampleEnd Ptr past the last sample node.
 *  @param tolerance Absolute tolerance to check ansatz function values against.
 */
template <class TScalarExpression, unsigned Dimension, concepts::CallableWith<BoundaryID,Size> TFunctor>
void scanConnectivities(Ref<const maths::AnsatzSpace<TScalarExpression,Dimension>> r_ansatzSpace,
                        TFunctor&& r_functor,
                        Ptr<const typename TScalarExpression::Value> p_sampleBegin,
                        Ptr<const typename TScalarExpression::Value> p_sampleEnd,
                        typename TScalarExpression::Value tolerance);



/** @brief Utility class for matching ansatz functions on different boundaries to preserve continuity.
 *  @details Example in 2D with ansatz functions generated from the outer product of the following set:
 *           @f[
 *              \begin{align}
 *                  \mathcal{N}_0(x) &= \frac{1 + x}{2}  \\
 *                  \mathcal{N}_1(x) &= \frac{1 - x}{2}  \\
 *                  \mathcal{N}_2(x) &= (1 + x) (1 - x)
 *              \end{align}
 *           @f]
 *           The resulting ansatz functions in 2D local space:
 *           @f[
 *              [N_i(\xi,\eta)] = \frac{1}{4}
 *              \begin{bmatrix}
 *                  (1+\xi) (1+\eta)         \\
 *                  (1-\xi) (1+\eta)         \\
 *                  2 (1-\xi^2) (1+\eta)     \\
 *                  (1+\xi) (1-\eta)         \\
 *                  (1-\xi) (1-\eta)         \\
 *                  2 (1-\xi^2) (1-\eta)     \\
 *                  2 (1+\xi) (1-\eta^2)     \\
 *                  2 (1-\xi) (1-\eta^2)     \\
 *                  4 (1-\xi^2) (1-\eta^2)
 *              \end{bmatrix}
 *           @f]
 *           @code
 *                     eta
 *
 *                      ^
 *                      |
 *               + ---- 1 ---- +
 *               |      |      |
 *               |      |      |
 *           -- -1 ---- + ---- 1 -- >   xi
 *               |      |      |
 *               |      |      |
 *               + --- -1 ---- +
 *                      |
 *                      |
 *           @endcode
 */
template <class TScalarExpression, unsigned Dimension>
class AnsatzMap : public Kernel<Dimension,typename TScalarExpression::Value>
{
private:
    using BaseTraits = Kernel<Dimension,typename TScalarExpression::Value>;

public:
    using typename BaseTraits::Value;

public:
    AnsatzMap(Ref<const maths::AnsatzSpace<TScalarExpression,Dimension>> r_ansatzSpace,
              std::span<const Value> samples,
              utils::Comparison<Value> comparison);

    template <class TOutputIt>
    void getPairs(BoundaryID boundary,
                  TOutputIt it_output) const noexcept;

    Size getPairCount(BoundaryID boundary) const noexcept;

private:
    tsl::robin_map<
        unsigned,                           // <== dimension index
        DynamicArray<std::pair<Size,Size>>  // <== list of coincident ansatz function indices on opposite boundaries
    > _connectivityMap;
}; // class AnsatzMap



template <class TScalarExpression, unsigned Dimension>
AnsatzMap<TScalarExpression,Dimension> makeAnsatzMap(Ref<const maths::AnsatzSpace<TScalarExpression,Dimension>> r_ansatzSpace,
                                                     std::span<const typename TScalarExpression::Value> samples,
                                                     utils::Comparison<typename TScalarExpression::Value> comparison);


} // namespace cie::fem

#include "packages/graph/impl/connectivity_impl.hpp"

#endif
