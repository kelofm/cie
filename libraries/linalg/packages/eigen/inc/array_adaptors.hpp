// --- External Includes ---
#include "Eigen/Dense"

// --- Utility Includes ---
#include "packages/compile_time/packages/concepts/inc/basic_concepts.hpp"
#include "packages/stl_extension/inc/StaticArray.hpp"
#include <packages/stl_extension/inc/DynamicArray.hpp>


namespace cie::linalg {


template <concepts::Numeric TValue, Size ArraySize>
Eigen::Map<Eigen::Matrix<TValue,ArraySize,1>> makeEigenAdaptor(Ref<StaticArray<TValue,ArraySize>> r_array);


template <concepts::Numeric TValue, Size ArraySize>
Eigen::Map<const Eigen::Matrix<TValue,ArraySize,1>> makeEigenAdaptor(Ref<const StaticArray<TValue,ArraySize>> r_array);


template <concepts::Numeric TValue, class ...TArgs>
Eigen::Map<Eigen::Matrix<TValue,Eigen::Dynamic,1>> makeEigenAdaptor(Ref<DynamicArray<TValue,TArgs...>> r_array);


template <concepts::Numeric TValue, class ...TArgs>
Eigen::Map<const Eigen::Matrix<TValue,Eigen::Dynamic,1>> makeEigenAdaptor(Ref<const DynamicArray<TValue,TArgs...>> r_array);


} // namespace cie::linalg
