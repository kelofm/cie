// --- FEM Includes ---
#include "packages/maths/inc/BSpline.hpp"

// --- GEO Includes ---
#include "packages/spline/inc/curve.hpp"

// --- STL Includes ---
#include <algorithm>
#include <format>


namespace cie::fem::maths {


template <class T, unsigned P, unsigned D>
BSpline<T,P,D>::BSpline(std::span<const std::span<const T>,PhysicalDimension> controlPointCoordinates,
                        std::span<const T> knots)
requires (ParametricDimension == 1u)
    : BSpline() {
    CIE_BEGIN_EXCEPTION_TRACING
    this->_controlPointCount.front() = controlPointCoordinates.front().size();
    _data.resize(this->controlPointCount(0u) * PhysicalDimension + knots.size());

    for (unsigned iDimension=0u; iDimension<PhysicalDimension; ++iDimension) {
        CIE_CHECK(
            controlPointCoordinates[iDimension].size() == this->controlPointCount(0),
            std::format(
                "inconsistent number of control point coordinates at dimension {}: expecting {} but got {}",
                iDimension,
                this->controlPointCount(0),
                controlPointCoordinates[iDimension].size()))
        std::copy_n(
            controlPointCoordinates[iDimension].data(),
            this->controlPointCount(0),
            this->controlPointCoordinates(iDimension).data());
    }

    if (!knots.empty()) {
        const T knotBegin = knots.front();
        const T knotEnd = knots.back();
        CIE_CHECK(knotEnd - knotBegin != static_cast<T>(0), "empty knot span")
        std::transform(
            knots.begin(),
            knots.end(),
            this->knots().data(),
            [knotBegin, knotEnd] (T knot) {
                return static_cast<T>(-1) + static_cast<T>(2) * (knot - knotBegin) / (knotEnd - knotBegin);
            });
    }

    CIE_END_EXCEPTION_TRACING
}


template <class T, unsigned P, unsigned D>
void BSpline<T,P,D>::evaluate(ConstSpan in, Span out) const {
    CIE_OUT_OF_RANGE_CHECK(ParametricDimension * out.size() == PhysicalDimension * in.size());
    if constexpr (ParametricDimension == 1u && PhysicalDimension == 2u) {
        geo::evaluateBSplineCurve<T,D>(
            in,
            {this->controlPointCoordinates(0), this->controlPointCoordinates(1)},
            this->knots(),
            {
                std::span<T>(
                    out.data(),
                    in.size()),
                std::span<T>(
                    out.data() + in.size(),
                    in.size())
            });
    } else {
        static_assert(
            std::is_same_v<T,void>,
            "spline not implemented for this dimension");
    }
}


template class BSpline<float, 1u, 2u>;
template class BSpline<double, 1u, 2u>;


} // namespace cie::fem::maths
