#ifndef CIE_LINALG_CONCEPTS_HPP
#define CIE_LINALG_CONCEPTS_HPP

// --- External Includes ---
#include <Eigen/Core>

// --- Utility Includes ---
#include "packages/compile_time/packages/concepts/inc/basic_concepts.hpp"


namespace cie::concepts {


namespace detail {

template <int I>
concept Dynamic
= (I == Eigen::Dynamic);

template <int I>
concept Static
= (!Dynamic<I>);

template <class T>
concept HasCIETraitTags
= requires()
{
    { T::RowTag };
    { T::ColumnTag };
};

template <class T>
concept HasEigenTraitTags
= requires()
{
    { T::RowsAtCompileTime };
    { T::ColsAtCompileTime };
};

template <class T>
struct MatrixTraitTags
{};

template <HasCIETraitTags T>
struct MatrixTraitTags<T>
{
    static const int RowTag = T::RowTag;
    static const int ColumnTag = T::ColumnTag;
};

template <HasEigenTraitTags T>
struct MatrixTraitTags<T>
{
    static const int RowTag    = T::RowsAtCompileTime;
    static const int ColumnTag = T::ColsAtCompileTime;
};

} // namespace detail


template <class T>
concept DynamicRowSize
= detail::Dynamic<detail::MatrixTraitTags<T>::RowTag>;


template <class T>
concept StaticRowSize
= (!DynamicRowSize<T>);


template <class T>
concept DynamicColumnSize
= detail::Dynamic<detail::MatrixTraitTags<T>::ColumnTag>;


template <class T>
concept StaticColumnSize
= (!DynamicColumnSize<T>);


template <class T>
concept DynamicMatrix
=  DynamicRowSize<T>
&& DynamicColumnSize<T>;


template <class T>
concept StaticMatrix
=  StaticRowSize<T>
&& StaticColumnSize<T>;


template <class T>
concept StaticArray
=  (StaticRowSize<T> && (detail::MatrixTraitTags<T>::ColumnTag == 1))
|| (StaticColumnSize<T> && (detail::MatrixTraitTags<T>::RowTag == 1));


template <class T>
concept DynamicArray
=  (DynamicRowSize<T> && (detail::MatrixTraitTags<T>::ColumnTag == 1))
|| (DynamicColumnSize<T> && (detail::MatrixTraitTags<T>::RowTag == 1));


} // namespace cie::concepts

#endif