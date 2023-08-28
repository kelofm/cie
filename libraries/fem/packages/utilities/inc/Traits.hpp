#ifndef CIE_FEM_UTILITIES_TRAITS_HPP
#define CIE_FEM_UTILITIES_TRAITS_HPP

// --- Linalg Includes ---
#include "packages/matrix/inc/concepts.hpp"

// --- Internal Includes ---
#include "packages/utilities/inc/kernel.hpp"


namespace cie::fem {


namespace detail {


template <class T>
concept HasTraits
= requires ()
{
    typename T::number_type;
    typename T::dynamic_array;
    typename T::Point;
    typename T::dense;
    typename T::sparse;
};

template <class T>
struct SizeTraits
{
    static const int RowTag = 0;
    static const int ColumnTag = 0;
};

template <class T>
requires (concepts::detail::HasCIETraitTags<T> || concepts::detail::HasEigenTraitTags<T>)
struct SizeTraits<T> : public concepts::detail::MatrixTraitTags<T>
{};


} // namespace detail


///@addtogroup fem
///@{


/// Invalid general Traits template
template <class ...Args>
struct Traits : public Kernel<0, void>
{};


template <Size Dimension, concepts::Numeric NT>
struct Traits<Dimension,NT> : public Kernel<Dimension,NT>
{};


template <detail::HasTraits T>
struct Traits<T> : public Kernel<T::Dimension, typename T::number_type>,
                   public detail::SizeTraits<T>
{};



///@}

} // namespace cie::fem


#endif