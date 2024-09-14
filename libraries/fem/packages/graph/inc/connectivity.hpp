#ifndef CIE_FEM_GRAPH_CONNECTIVITY_HPP
#define CIE_FEM_GRAPH_CONNECTIVITY_HPP

// --- External Includes ---
#include "tsl/robin_map.h"

// --- FEM Includes ---
#include "packages/maths/inc/Expression.hpp"
#include "packages/graph/inc/BoundaryID.hpp"
#include "packages/utilities/inc/kernel.hpp"

// --- Utility Includes ---
#include "packages/compile_time/packages/concepts/inc/functional.hpp"
#include "packages/compile_time/packages/concepts/inc/iterator_concepts.hpp"
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
 *  @param rAnsatzSpace @ref maths::AnsatzSpace to scan the functions of.
 *  @param rFunctor Functor that gets called with each @ref BoundaryID and non-vanishing
 *                   ansatz function index.
 *  @param pSampleBegin Ptr to the beginning of the array of sample nodes to evaluate
 *                       the ansatz functions at.
 *  @param pSampleEnd Ptr past the last sample node.
 *  @param tolerance Absolute tolerance to check ansatz function values against.
 */
template <maths::Expression TAnsatzSpace, concepts::CallableWith<BoundaryID,Size> TFunctor>
void scanConnectivities(Ref<const TAnsatzSpace> rAnsatzSpace,
                        TFunctor&& rFunctor,
                        Ptr<const typename TAnsatzSpace::Value> pSampleBegin,
                        Ptr<const typename TAnsatzSpace::Value> pSampleEnd,
                        typename TAnsatzSpace::Value tolerance);



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
 *              [NI(\xi,\eta)] = \frac{1}{4}
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
template <class TValue>
class AnsatzMap
{
public:
    template <maths::Expression TAnsatzSpace>
    requires (std::is_same_v<typename TAnsatzSpace::Value,TValue>)
    AnsatzMap(Ref<const TAnsatzSpace> rAnsatzSpace,
              std::span<const TValue> samples,
              utils::Comparison<TValue> comparison);

    template <concepts::OutputIterator<std::pair<Size,Size>> TOutputIt>
    void getPairs(BoundaryID first,
                  BoundaryID second,
                  TOutputIt itOutput) const noexcept;

    Size getPairCount(BoundaryID boundary) const noexcept;

private:
    tsl::robin_map<
        std::pair<BoundaryID,BoundaryID>,   // <== dimension index
        DynamicArray<std::pair<Size,Size>>  // <== list of coincident ansatz function indices on opposite boundaries
    > _connectivityMap;
}; // class AnsatzMap



template <maths::Expression TAnsatzSpace>
AnsatzMap<typename TAnsatzSpace::Value> makeAnsatzMap(Ref<const TAnsatzSpace> rAnsatzSpace,
                                                      std::span<const typename TAnsatzSpace::Value> samples,
                                                      utils::Comparison<typename TAnsatzSpace::Value> comparison);


} // namespace cie::fem

#include "packages/graph/impl/connectivity_impl.hpp"

#endif
