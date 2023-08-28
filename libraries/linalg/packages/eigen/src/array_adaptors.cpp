// --- LinAlg Includes ---
#include "packages/eigen/inc/array_adaptors.hpp"


namespace cie::linalg {


template <concepts::Numeric TValue, Size ArraySize>
Eigen::Map<Eigen::Matrix<TValue,ArraySize,1>> makeEigenAdaptor(Ref<StaticArray<TValue,ArraySize>> r_array)
{
    return Eigen::Map<Eigen::Matrix<TValue,ArraySize,1>>(r_array.data());
}


template <concepts::Numeric TValue, Size ArraySize>
Eigen::Map<const Eigen::Matrix<TValue,ArraySize,1>> makeEigenAdaptor(Ref<const StaticArray<TValue,ArraySize>> r_array)
{
    return Eigen::Map<const Eigen::Matrix<TValue,ArraySize,1>>(r_array.data());
}


template <concepts::Numeric TValue, class ...TArgs>
Eigen::Map<Eigen::Matrix<TValue,Eigen::Dynamic,1>> makeEigenAdaptor(Ref<DynamicArray<TValue,TArgs...>> r_array)
{
    return Eigen::Map<Eigen::Matrix<TValue,Eigen::Dynamic,1>>(
        r_array.data(),
        r_array.size()
    );
}


template <concepts::Numeric TValue, class ...TArgs>
Eigen::Map<const Eigen::Matrix<TValue,Eigen::Dynamic,1>> makeEigenAdaptor(Ref<const DynamicArray<TValue,TArgs...>> r_array)
{
    return Eigen::Map<const Eigen::Matrix<TValue,Eigen::Dynamic,1>>(
        r_array.data(),
        r_array.size()
    );
}


#define CIE_LINALG_INSTANTIATE_STATIC_EIGEN_ADAPTORS(TValue, ArraySize)             \
    template <> Eigen::Map<Eigen::Matrix<TValue,ArraySize,1>>                       \
    makeEigenAdaptor<TValue,ArraySize>(Ref<StaticArray<TValue,ArraySize>>);         \
    template <> Eigen::Map<const Eigen::Matrix<TValue,ArraySize,1>>                 \
    makeEigenAdaptor<TValue,ArraySize>(Ref<const StaticArray<TValue,ArraySize>>)


#define CIE_LINALG_INSTANTIATE_EIGEN_ADAPTORS(TValue)                               \
    CIE_LINALG_INSTANTIATE_STATIC_EIGEN_ADAPTORS(TValue, 1);                        \
    CIE_LINALG_INSTANTIATE_STATIC_EIGEN_ADAPTORS(TValue, 2);                        \
    CIE_LINALG_INSTANTIATE_STATIC_EIGEN_ADAPTORS(TValue, 3);                        \
    CIE_LINALG_INSTANTIATE_STATIC_EIGEN_ADAPTORS(TValue, 4);                        \
    CIE_LINALG_INSTANTIATE_STATIC_EIGEN_ADAPTORS(TValue, 6);                        \
    CIE_LINALG_INSTANTIATE_STATIC_EIGEN_ADAPTORS(TValue, 9);                        \
    template <> Eigen::Map<Eigen::Matrix<TValue,Eigen::Dynamic,1>>                  \
    makeEigenAdaptor(Ref<DynamicArray<TValue>>);                                    \
    template <> Eigen::Map<const Eigen::Matrix<TValue,Eigen::Dynamic,1>>            \
    makeEigenAdaptor(Ref<const DynamicArray<TValue>>)


CIE_LINALG_INSTANTIATE_EIGEN_ADAPTORS(unsigned);

CIE_LINALG_INSTANTIATE_EIGEN_ADAPTORS(int);

CIE_LINALG_INSTANTIATE_EIGEN_ADAPTORS(float);

CIE_LINALG_INSTANTIATE_EIGEN_ADAPTORS(double);


#undef CIE_LINALG_INSTANTIATE_EIGEN_ADAPTORS

#undef CIE_LINALG_INSTANTIATE_STATIC_EIGEN_ADAPTORS


} // namespace cie::linalg
