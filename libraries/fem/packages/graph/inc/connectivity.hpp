#ifndef CIE_FEM_GRAPH_CONNECTIVITY_HPP
#define CIE_FEM_GRAPH_CONNECTIVITY_HPP

// --- FEM Includes ---
#include "packages/maths/inc/AnsatzSpace.hpp"
#include "packages/graph/inc/BoundaryID.hpp"

// --- Utility Includes ---
#include "packages/compile_time/packages/concepts/inc/functional.hpp"


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
template <class TScalarExpression, unsigned Dimension, concepts::CallableWith<BoundaryID,unsigned> TFunctor>
void scanConnectivities(Ref<const maths::AnsatzSpace<TScalarExpression,Dimension>> r_ansatzSpace,
                        TFunctor&& r_functor,
                        Ptr<const typename TScalarExpression::Value> p_sampleBegin,
                        Ptr<const typename TScalarExpression::Value> p_sampleEnd,
                        typename TScalarExpression::Value tolerance);


} // namespace cie::fem

#include "packages/graph/impl/connectivity_impl.hpp"

#endif
