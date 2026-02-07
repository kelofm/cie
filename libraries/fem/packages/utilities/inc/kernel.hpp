#ifndef CIE_FEM_UTILITIES_KERNEL_HPP
#define CIE_FEM_UTILITIES_KERNEL_HPP

// --- External Includes ---
#include "Eigen/Dense"

// --- Utility Includes ---
#include "packages/compile_time/packages/concepts/inc/basic_concepts.hpp"
#include "packages/compile_time/packages/parameter_pack/inc/Match.hpp"
#include "packages/types/inc/types.hpp"
#include "packages/stl_extension/inc/StrongTypeDef.hpp"
#include "packages/stl_extension/inc/StrongTypeDef_overloads.hpp"
#include "packages/stl_extension/inc/StaticArray.hpp"

// --- LinAlg Includes ---
#include "packages/matrix/inc/StaticEigenMatrix.hpp"
#include "packages/matrix/inc/DynamicEigenMatrix.hpp"
#include "packages/matrix/inc/SparseEigenMatrix.hpp"

// --- STL Includes ---
#include <span>


namespace cie::fem {


///@addtogroup fem
///@{

namespace detail {

struct ParametricSpaceTag {};

struct PhysicalSpaceTag {};

template <unsigned Dimension, concepts::Numeric NT>
using Point = linalg::StaticEigenMatrix<NT,Dimension,1>;

template <unsigned Dimension, concepts::Numeric NT>
using LocalPoint = utils::StrongTypeDef<Point<Dimension,NT>,ParametricSpaceTag>;

template <unsigned Dimension, concepts::Numeric NT>
using GlobalPoint = utils::StrongTypeDef<Point<Dimension,NT>,PhysicalSpaceTag>;

} // namespace detail


template <concepts::Numeric T>
using ParametricCoordinate = utils::StrongTypeDef<T,detail::ParametricSpaceTag>;

template <concepts::Numeric T>
using PhysicalCoordinate = utils::StrongTypeDef<T,detail::PhysicalSpaceTag>;


template <unsigned Dimension, concepts::Numeric NT>
struct Kernel
{
    static const unsigned dimension = Dimension;
    using Value                     = NT;

    using number_type               = NT;
    using dynamic_array             = linalg::EigenMatrix<Eigen::Matrix<NT,Eigen::Dynamic,1>>;

    template <Size ArraySize>
    using static_array              = linalg::StaticEigenMatrix<NT,ArraySize,1>;

    using ParametricCoordinate      = ParametricCoordinate<NT>;
    using PhysicalCoordinate        = PhysicalCoordinate<NT>;

    using Point                     = detail::Point<Dimension,NT>;
    using LocalPoint                = detail::LocalPoint<Dimension,NT>;
    using GlobalPoint               = detail::GlobalPoint<Dimension,NT>;

    struct dense {
        using dynamic_matrix = linalg::DynamicEigenMatrix<NT>;

        using DynamicAdaptor = Eigen::Map<dynamic_matrix>;

        template <Size RowSize, Size ColumnSize>
        using static_matrix = linalg::StaticEigenMatrix<NT,RowSize,ColumnSize>;
    };

    struct sparse {
        using dynamic_matrix = linalg::SparseEigenMatrix<NT>;

        template <Size RowSize, Size ColumnSize>
        using static_matrix = void; // dummy
    };

    /// @brief Strip local system type information from a view over coordinates.
    /// @param rSpan Input view over local coordinates.
    template <std::size_t SpanSize>
    static constexpr std::span<const Value,SpanSize> decay(const std::span<const ParametricCoordinate,SpanSize>& rSpan) noexcept {
        return std::span<const Value,SpanSize>(
            reinterpret_cast<const Value*>(rSpan.data()),
            SpanSize);
    }

    /// @brief Strip local system type information from a view over coordinates.
    /// @param rSpan Input view over local coordinates.
    template <std::size_t SpanSize>
    static constexpr std::span<Value,SpanSize> decay(const std::span<ParametricCoordinate,SpanSize>& rSpan) noexcept {
        return std::span<Value,SpanSize>(
            reinterpret_cast<Value*>(rSpan.data()),
            SpanSize);
    }

    /// @brief Strip global system type information from a view over coordinates.
    /// @param rSpan Input view over global coordinates.
    template <std::size_t SpanSize>
    static constexpr std::span<const Value,SpanSize> decay(const std::span<const PhysicalCoordinate,SpanSize>& rSpan) noexcept {
        return std::span<const Value,SpanSize>(
            reinterpret_cast<const Value*>(rSpan.data()),
            SpanSize);
    }

    /// @brief Strip global system type information from a view over coordinates.
    /// @param rSpan Input view over global coordinates.
    template <std::size_t SpanSize>
    static constexpr std::span<Value,SpanSize> decay(const std::span<PhysicalCoordinate,SpanSize>& rSpan) noexcept {
        return std::span<Value,SpanSize>(
            reinterpret_cast<Value*>(rSpan.data()),
            SpanSize);
    }

    template <class T, std::size_t SpanSize>
    requires ct::Match<T>::template Any<ParametricCoordinate,PhysicalCoordinate>
    static constexpr std::span<const T,SpanSize> cast(const std::span<const Value,SpanSize>& rSpan) noexcept {
        return std::span<const T,SpanSize>(
            reinterpret_cast<const T*>(rSpan.data()),
            SpanSize);
    }

    template <class T, std::size_t SpanSize>
    requires ct::Match<T>::template Any<ParametricCoordinate,PhysicalCoordinate>
    static constexpr std::span<T,SpanSize> cast(const std::span<Value,SpanSize>& rSpan) noexcept {
        return std::span<T,SpanSize>(
            reinterpret_cast<T*>(rSpan.data()),
            SpanSize);
    }

    template <class TIn, std::size_t ArraySize>
    requires (ct::Match<TIn>::template Any<Value,ParametricCoordinate,PhysicalCoordinate>)
    static constexpr std::span<const TIn,ArraySize> view(const StaticArray<TIn,ArraySize>& rArray) noexcept {
        return std::span<const TIn,ArraySize>(
            reinterpret_cast<const TIn*>(rArray.data()),
            ArraySize);
    }

    template <class TIn, std::size_t ArraySize>
    requires ct::Match<TIn>::template Any<Value,ParametricCoordinate,PhysicalCoordinate>
    static constexpr std::span<TIn,ArraySize> view(StaticArray<TIn,ArraySize>& rArray) noexcept {
        return std::span<TIn,ArraySize>(
            reinterpret_cast<TIn*>(rArray.data()),
            ArraySize);
    }

    template <class TIn, std::size_t ArraySize>
    requires ct::Match<TIn>::template Any<Value,ParametricCoordinate,PhysicalCoordinate>
    static constexpr std::span<const Value,ArraySize> decayView(const StaticArray<TIn,ArraySize>& rArray) noexcept {
        return Kernel::decay(Kernel::view(rArray));
    }

    template <class TIn, std::size_t ArraySize>
    requires ct::Match<TIn>::template Any<Value,ParametricCoordinate,PhysicalCoordinate>
    static constexpr std::span<Value,ArraySize> decayView(StaticArray<TIn,ArraySize>& rArray) noexcept {
        return Kernel::decay(Kernel::view(rArray));
    }

    template <class TOut, class TIn, std::size_t ArraySize>
    requires ct::Match<TIn,TOut>::template All<Value,ParametricCoordinate,PhysicalCoordinate>
    static constexpr std::span<const TOut,ArraySize> castView(const StaticArray<TIn,ArraySize>& rArray) noexcept {
        return Kernel::cast<TOut>(Kernel::view(rArray));
    }

    template <class TOut, class TIn, std::size_t ArraySize>
    requires ct::Match<TIn,TOut>::template All<Value,ParametricCoordinate,PhysicalCoordinate>
    static constexpr std::span<TOut,ArraySize> castView(StaticArray<TIn,ArraySize>& rArray) noexcept {
        return Kernel::cast<TOut>(Kernel::view(rArray));
    }

    template <class TIn, std::size_t ArraySize>
    requires (ct::Match<TIn>::template Any<Value,ParametricCoordinate,PhysicalCoordinate>)
    static constexpr std::span<const TIn,ArraySize> view(const std::array<TIn,ArraySize>& rArray) noexcept {
        return std::span<const TIn,ArraySize>(
            reinterpret_cast<const TIn*>(rArray.data()),
            ArraySize);
    }

    template <class TIn, std::size_t ArraySize>
    requires ct::Match<TIn>::template Any<Value,ParametricCoordinate,PhysicalCoordinate>
    static constexpr std::span<TIn,ArraySize> view(std::array<TIn,ArraySize>& rArray) noexcept {
        return std::span<TIn,ArraySize>(
            reinterpret_cast<TIn*>(rArray.data()),
            ArraySize);
    }

    template <class TIn, std::size_t ArraySize>
    requires ct::Match<TIn>::template Any<Value,ParametricCoordinate,PhysicalCoordinate>
    static constexpr std::span<const Value,ArraySize> decayView(const std::array<TIn,ArraySize>& rArray) noexcept {
        return Kernel::decay(Kernel::view(rArray));
    }

    template <class TIn, std::size_t ArraySize>
    requires ct::Match<TIn>::template Any<Value,ParametricCoordinate,PhysicalCoordinate>
    static constexpr std::span<Value,ArraySize> decayView(std::array<TIn,ArraySize>& rArray) noexcept {
        return Kernel::decay(Kernel::view(rArray));
    }

    template <class TOut, class TIn, std::size_t ArraySize>
    requires ct::Match<TIn,TOut>::template All<Value,ParametricCoordinate,PhysicalCoordinate>
    static constexpr std::span<const TOut,ArraySize> castView(const std::array<TIn,ArraySize>& rArray) noexcept {
        return Kernel::cast<TOut>(Kernel::view(rArray));
    }

    template <class TOut, class TIn, std::size_t ArraySize>
    requires ct::Match<TIn,TOut>::template All<Value,ParametricCoordinate,PhysicalCoordinate>
    static constexpr std::span<TOut,ArraySize> castView(std::array<TIn,ArraySize>& rArray) noexcept {
        return Kernel::cast<TOut>(Kernel::view(rArray));
    }
}; // class Kernel

///@}

} // namespace cie::fem


#endif