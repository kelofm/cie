#pragma once

// --- FEM Includes ---
#include "packages/numeric/inc/QuadraturePoint.hpp"
#include "packages/numeric/inc/Cell.hpp"
#include "packages/maths/inc/IdentityTransform.hpp"
#include "packages/graph/inc/Graph.hpp"

// --- Utility Includes ---
#include "packages/compile_time/packages/concepts/inc/basic_concepts.hpp"

// --- STL Includes ---
#include <type_traits>


namespace cie::fem {


template <class T, class TQuadraturePointData = void>
concept QuadraturePointFactoryLike
=  cie::concepts::UnsignedInteger<std::remove_const_t<decltype(T::Dimension)>>
&& cie::concepts::Numeric<typename T::Value>
&& requires(T& rInstance, std::span<QuadraturePoint<T::Dimension,typename T::Value,TQuadraturePointData>> quadraturePoints) {
    {rInstance(quadraturePoints)} -> cie::concepts::UnsignedInteger;
}; // concept QuadraturePointFactoryLike


template <
    class T,
    class TMesh = Graph<CellBase<1,float,maths::IdentityTransform<float,1>>,void,void>,
    class TCell = typename TMesh::Vertex::Data,
    class TQuadraturePointData = void>
concept QuadratureRuleFactoryLike
=  GraphLike<TMesh>
&& CellLike<TCell>
&& requires (const T& constInstance, const TMesh& rMesh, const TCell& rCell) {
    {constInstance(rMesh, rCell)} -> QuadraturePointFactoryLike<TQuadraturePointData>;
}; // concept QuadratureRuleFactoryLike



template <
    unsigned Dim,
    concepts::Numeric TValue,
    class TQuadraturePointData = void>
class OuterProductQuadraturePointFactory {
public:
    static constexpr inline unsigned Dimension = Dim;

    using Value = TValue;

    constexpr OuterProductQuadraturePointFactory() noexcept;

    constexpr explicit OuterProductQuadraturePointFactory(std::span<const QuadraturePoint<1u,TValue,TQuadraturePointData>> base) noexcept;

    unsigned operator()(std::span<QuadraturePoint<Dimension,TValue,TQuadraturePointData>> out) noexcept;

private:
    bool _done;

    std::span<const QuadraturePoint<1u,TValue,TQuadraturePointData>> _base;

    StaticArray<unsigned,Dimension> _state;
}; // class OuterProductQuadraturePointFactory


template <
    unsigned Dim,
    concepts::Numeric TValue,
    class TQuadraturePointData = void>
class CachedQuadraturePointFactory {
public:
    static constexpr inline unsigned Dimension = Dim;

    using Value = TValue;

    CachedQuadraturePointFactory() noexcept;

    explicit CachedQuadraturePointFactory(std::span<const QuadraturePoint<Dimension,Value,TQuadraturePointData>> cache) noexcept;

    unsigned operator()(std::span<QuadraturePoint<Dimension,TValue>> out) noexcept;

private:
    std::span<const QuadraturePoint<Dimension,Value,TQuadraturePointData>> _cache;
};


} // namespace cie::fem

#include "packages/numeric/impl/QuadraturePointFactory_impl.hpp"
