#pragma once

// --- External Includes ---
#include "tsl/robin_map.h"

// --- FEM Includes ---
#include "packages/maths/inc/Expression.hpp"
#include "packages/graph/inc/BoundaryID.hpp"
#include "packages/graph/inc/OrientedBoundary.hpp"
#include "packages/utilities/inc/kernel.hpp"

// --- Utility Includes ---
#include "packages/compile_time/packages/concepts/inc/functional.hpp"
#include "packages/compile_time/packages/concepts/inc/iterator_concepts.hpp"
#include "packages/stl_extension/inc/DynamicArray.hpp"
#include "packages/maths/inc/Comparison.hpp"
#include "packages/stl_extension/inc/Hash.hpp"

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
template <maths::Expression TAnsatzSpace, cie::concepts::CallableWith<BoundaryID,Size> TFunctor>
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
template <unsigned Dim>
class AnsatzMap {
private:
    struct SymmetricHash {
        auto operator()(const std::pair<OrientedBoundary<Dim>,OrientedBoundary<Dim>>& rPair) const noexcept {
            utils::Hash<std::pair<OrientedBoundary<Dim>,OrientedBoundary<Dim>>> hasher;
            if (rPair.second < rPair.first) return hasher(std::make_pair(rPair.second, rPair.first));
            else return hasher(rPair);
        }
    }; // struct SymmetricHash

    struct SymmetricEquality {
        auto operator()(const std::pair<OrientedBoundary<Dim>,OrientedBoundary<Dim>>& rLhs,
                        const std::pair<OrientedBoundary<Dim>,OrientedBoundary<Dim>>& rRhs) const noexcept {
            return (rLhs.first == rRhs.first && rLhs.second == rRhs.second)
                || (rLhs.second == rRhs.first && rLhs.first == rRhs.second);
        }
    }; // struct SymmetricOrdering

    using ConnectivityMap = tsl::robin_map<
        std::pair<OrientedBoundary<Dim>,OrientedBoundary<Dim>>, ///< oriented boundary pair
        DynamicArray<std::pair<std::size_t,std::size_t>>,       ///< list of coincident ansatz function indices on the pair of boundaries
        SymmetricHash,
        SymmetricEquality
    >;

    struct AnsatzPairs {
        typename ConnectivityMap::const_iterator it;
        bool swap;
    }; // class AnsatzPairs

public:
    static constexpr inline unsigned Dimension = Dim;

    AnsatzMap() noexcept = default;

    template <maths::Expression TAnsatzSpace>
    AnsatzMap(
        Ref<const TAnsatzSpace> rAnsatzSpace,
        std::size_t integrationOrder,
        utils::Comparison<typename TAnsatzSpace::Value> comparison = utils::Comparison<typename TAnsatzSpace::Value>());

    /// @brief Find coincident ansatz functions in a direction.
    /// @see @ref pairCount
    /// @see @ref getPairs
    [[nodiscard]] AnsatzPairs findPairs(
        const OrientedBoundary<Dim> first,
        const OrientedBoundary<Dim> second) const noexcept;

    /// @brief Fetch the number of coincident ansatz functions in a direction.
    /// @see @ref findPairs
    /// @see @ref getPairs
    [[nodiscard]] std::size_t pairCount(AnsatzPairs pairs) const noexcept;

    /// @brief Fetch the indices of coincident ansatz functions in a direction.
    /// @see @ref findPairs
    /// @see @ref pairCount
    void getPairs(
        AnsatzPairs pairs,
        std::span<std::pair<std::size_t,std::size_t>> output) const;

    /// @brief Get the number of ansatz function in the @ref AnsatzSpace.
    [[nodiscard]] std::size_t ansatzCount() const noexcept;

private:
    std::size_t _ansatzCount;

    ConnectivityMap _topology;
}; // class AnsatzMap



template <maths::Expression TAnsatzSpace, class TValue>
requires std::is_same_v<TValue,typename TAnsatzSpace::Value>
AnsatzMap<TAnsatzSpace::Dimension>
makeAnsatzMap(
    Ref<const TAnsatzSpace> rAnsatzSpace,
    std::size_t integrationOrder,
    utils::Comparison<TValue> comparison = utils::Comparison<typename TAnsatzSpace::Value>());


} // namespace cie::fem

#include "packages/graph/impl/connectivity_impl.hpp"
