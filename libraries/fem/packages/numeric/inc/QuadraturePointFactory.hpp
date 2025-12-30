#pragma once

// --- FEM Includes ---
#include "packages/numeric/inc/QuadraturePoint.hpp"
#include "packages/numeric/inc/Cell.hpp"
#include "packages/maths/inc/IdentityTransform.hpp"

// --- Utility Includes ---
#include "packages/compile_time/packages/concepts/inc/basic_concepts.hpp"
#include <type_traits>


namespace cie::fem {


template <class T, class TCell = CellBase<1,float,maths::IdentityTransform<float,1>>>
concept QuadraturePointFactoryLike =
   CellLike<TCell>
&& std::is_same_v<std::remove_cvref_t<decltype(T::Dimension)>,unsigned>
&& ::cie::concepts::Numeric<typename T::Value>
&& requires(T& rInstance,
            const TCell& rCell,
            std::span<QuadraturePoint<T::Dimension,typename T::Value>> quadraturePoints) {
    {rInstance.generate(rCell, quadraturePoints)} -> cie::concepts::UnsignedInteger;
}; // concept QuadraturePointFactoryLike


template <unsigned Dim, concepts::Numeric TValue>
class OuterProductQuadraturePointFactory {
public:
    static constexpr inline unsigned Dimension = Dim;

    using Value = TValue;

    constexpr OuterProductQuadraturePointFactory() noexcept;

    constexpr explicit OuterProductQuadraturePointFactory(std::span<const QuadraturePoint<1u,TValue>> base) noexcept;

    template <CellLike TCell>
    unsigned generate(Ref<const TCell>,
                      std::span<QuadraturePoint<Dimension,TValue>> out) noexcept;

private:
    bool _done;

    std::span<const QuadraturePoint<1u,TValue>> _base;

    StaticArray<unsigned,Dimension> _state;
}; // class OuterProductQuadraturePointFactory


template <unsigned Dim, concepts::Numeric TValue>
class CachedQuadraturePointFactory {
public:
    static constexpr inline unsigned Dimension = Dimension;

    using Value = TValue;

    CachedQuadraturePointFactory() noexcept;

    CachedQuadraturePointFactory(std::span<const QuadraturePoint<Dimension,Value>> cache) noexcept;

    template <CellLike TCell>
    unsigned generate(Ref<const TCell>,
                      std::span<QuadraturePoint<Dimension,TValue>> out) noexcept;

private:
    std::span<const QuadraturePoint<Dimension,Value>> _cache;
};


} // namespace cie::fem

#include "packages/numeric/impl/QuadraturePointFactory_impl.hpp"
