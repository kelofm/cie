#ifndef CIE_FEM_UTILITIES_KERNEL_HPP
#define CIE_FEM_UTILITIES_KERNEL_HPP

// --- External Includes ---
#include "Eigen/Dense"

// --- Utility Includes ---
#include "packages/compile_time/packages/concepts/inc/basic_concepts.hpp"
#include "packages/types/inc/types.hpp"
#include "packages/stl_extension/inc/StaticArray.hpp"
#include "packages/stl_extension/inc/StrongTypeDef.hpp"

// --- LinAlg Includes ---
#include "packages/matrix/inc/StaticEigenMatrix.hpp"
#include "packages/matrix/inc/DynamicEigenMatrix.hpp"
#include "packages/matrix/inc/SparseEigenMatrix.hpp"

// --- STL Includes ---
#include <array>
#include <complex>


namespace cie::fem {


///@addtogroup fem
///@{

namespace detail {

template <Size Dimension, concepts::Numeric NT>
using Point = linalg::StaticEigenMatrix<NT,Dimension,1>;

struct LocalPointTag {};

template <Size Dimension, concepts::Numeric NT>
using LocalPoint = utils::StrongTypeDef<Point<Dimension,NT>,LocalPointTag>;

struct GlobalPointTag {};

template <Size Dimension, concepts::Numeric NT>
using GlobalPoint = utils::StrongTypeDef<Point<Dimension,NT>,GlobalPointTag>;

} // namespace detail


template <Size Dimension, concepts::Numeric NT>
struct Kernel
{
    static const Size dimension = Dimension;
    using number_type           = NT;
    using dynamic_array         = linalg::EigenMatrix<Eigen::Matrix<NT,Eigen::Dynamic,1>>;

    template <Size ArraySize>
    using static_array          = linalg::StaticEigenMatrix<NT,ArraySize,1>;

    using Point                 = detail::Point<Dimension,NT>;
    using LocalPoint            = detail::LocalPoint<Dimension,NT>;
    using GlobalPoint           = detail::GlobalPoint<Dimension,NT>;

    struct dense
    {
        using dynamic_matrix = linalg::DynamicEigenMatrix<NT>;

        template <Size RowSize, Size ColumnSize>
        using static_matrix = linalg::StaticEigenMatrix<NT,RowSize,ColumnSize>;
    };

    struct sparse
    {
        using dynamic_matrix = linalg::SparseEigenMatrix<NT>;

        template <Size RowSize, Size ColumnSize>
        using static_matrix = void; // dummy
    };
}; // class Kernel

///@}

} // namespace cie::fem


#endif