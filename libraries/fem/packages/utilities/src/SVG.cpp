// --- FEM Includes ---
#include "packages/utilities/inc/SVG.hpp"

// --- GEO Includes ---
#include "packages/spline/inc/interpolation.hpp"

// --- STL Includes ---
#include <format>


namespace cie::fem {


template <concepts::Numeric T>
SVG<T>::SVG(RightRef<std::vector<CurveVariant>> rCurves)
    : _curves(std::move(rCurves))
{}


template <concepts::Numeric T>
void SVG<T>::add(RightRef<CurveVariant> rCurve) {
    CIE_CHECK(
        !std::holds_alternative<std::monostate>(rCurve),
        "attempting to add uninitialized curve to SVG")
    _curves.emplace_back(std::move(rCurve));
}


template <concepts::Numeric T>
std::vector<std::array<T,2>>
SVG<T>::tesselate(std::size_t segmentsPerCurve) const {
    std::vector<std::array<T,2>> output;

    if (segmentsPerCurve) {
        CIE_BEGIN_EXCEPTION_TRACING
        output.reserve(_curves.size() * (segmentsPerCurve + 1));

        std::vector<T> in, out;
        in.resize(segmentsPerCurve + 1);
        out.resize(in.size() * 2);

        {
            const T d = static_cast<T>(2) / segmentsPerCurve;
            for (std::size_t iSample=0ul; iSample<in.size(); ++iSample)
                in[iSample] = static_cast<T>(-1) + iSample * d;
        }

        for (std::size_t iCurve=0ul; iCurve<_curves.size(); ++iCurve) {
            std::visit(
                [&in, &out] (const auto& rCurve) {
                    using TCurve = std::remove_cvref_t<decltype(rCurve)>;
                    CIE_BEGIN_EXCEPTION_TRACING
                    if constexpr (std::is_same_v<TCurve,maths::BSpline<T,1,2>>) {
                        rCurve.evaluate(in, out);
                    } else if constexpr (std::is_same_v<TCurve,maths::AffineEmbedding<T,1,2>>) {
                        for (std::size_t iIn=0ul; iIn<in.size(); ++iIn) {
                            std::array<T,2> buffer;
                            rCurve.evaluate(
                                {in.data() + iIn, 1},
                                buffer);
                            out[iIn] = buffer[0];
                            out[iIn + in.size()] = buffer[1];
                        }
                    } else {
                        CIE_THROW(Exception, "uninitialized curve segment")
                    }
                    CIE_END_EXCEPTION_TRACING
                }, _curves[iCurve]);
            std::transform(
                out.begin(),
                out.begin() + in.size(),
                out.begin() + in.size(),
                std::back_inserter(output),
                [] (T x, T y) {return std::array<T,2> {x, y};});
        } // for rCurve in _curves
        CIE_END_EXCEPTION_TRACING
    } // if segmentsPerCurve

    return output;
}


template class SVG<float>;

template class SVG<double>;


} // namespace cie::fem


namespace cie::fem::io {


template <class T>
GraphML::Deserializer<fem::SVG<T>>::Deserializer(
        Ref<fem::SVG<T>> rInstance,
        Ref<SAXHandler> rSAX,
        Ref<const maths::AffineTransform<T,2>> rTransform,
        std::string_view parseMode)
    : GraphML::DeserializerBase<fem::SVG<T>>(rInstance, rSAX),
      _mode(parseMode),
      _attributes(),
      _transform(rTransform)
{}


template <class T>
void GraphML::Deserializer<fem::SVG<T>>::onElementBegin(
        Ptr<void> pThis,
        std::string_view elementName,
        std::span<GraphML::AttributePair> attributes) {
    Ref<Deserializer> rThis = *static_cast<Ptr<Deserializer>>(pThis);
    bool isExpectedElement = false;
    isExpectedElement |= elementName == "svg";
    isExpectedElement |= elementName == "path";
    isExpectedElement |= elementName == "defs";
    CIE_CHECK(
        isExpectedElement,
        std::format(
            "unexpected element '{}'",
            elementName))

    if (rThis._mode == "svg") {
        if (elementName == "svg") {
            std::array<typename Kernel<2,T>::Point,3> transformed;
            transformed[0][0] = static_cast<T>(-1);
            transformed[0][1] = static_cast<T>( 1);
            transformed[1][0] = static_cast<T>( 1);
            transformed[1][1] = static_cast<T>( 1);
            transformed[2][0] = static_cast<T>(-1);
            transformed[2][1] = static_cast<T>(-1);
            rThis._transform = maths::AffineTransform<T,2>(transformed);
        } else {
            auto pDeserializer = Deserializer::make(
                rThis.instance(),
                rThis.sax(),
                elementName,
                rThis._transform,
                elementName);
            pDeserializer->onElementBegin(
                pDeserializer,
                elementName,
                attributes);
            rThis.sax().push({
                pDeserializer,
                Deserializer::onElementBegin,
                Deserializer::onText,
                Deserializer::onElementEnd});
        }
    } else {
        for (const auto& [attributeName, attributeValue] : attributes) {
            Ref<std::string> rAttributeValue = rThis._attributes.emplace(attributeName, std::string()).first->second;
            rAttributeValue.insert(
                rAttributeValue.end(),
                attributeValue.begin(),
                attributeValue.end());
        } // for attributeName, attributeValue in attributes
    }
}


template <class T>
void GraphML::Deserializer<fem::SVG<T>>::onText(
        Ptr<void>,
        std::string_view t) {
    if (std::find_if(t.begin(), t.end(), [](auto c) {return !std::isspace(c);}) != t.end())
        CIE_THROW(
            Exception,
            std::format("Unexpected text data on a <graphml> element: '{}'.", t))
}


template <class T, unsigned D>
maths::BSpline<T,1,D> makeBSpline(std::span<const std::vector<T>,D> controlPoints) {
    CIE_CHECK(
        controlPoints.size() == 2,
        std::format(
            "unexpected number of dimensions {} for a 2D BSpline",
            controlPoints.size()))
    CIE_CHECK(
        4 <= controlPoints.front().size(),
        std::format(
            "insufficient number of control points ({})",
            controlPoints.front().size()))

    if (4 == controlPoints.front().size()) {
        CIE_BEGIN_EXCEPTION_TRACING
        std::array<std::span<const T>,2> spans {
            controlPoints[0],
            controlPoints[1]};
        std::array<T,8> knots;
        for (std::size_t iKnot=0ul; iKnot<4ul; ++iKnot)
            knots[iKnot] = static_cast<T>(-1);
        for (std::size_t iKnot=4ul; iKnot<8ul; ++iKnot)
            knots[iKnot] = static_cast<T>( 1);
        return maths::BSpline<T,1,2>(spans, knots);
        CIE_END_EXCEPTION_TRACING
    } else {
        CIE_BEGIN_EXCEPTION_TRACING
        // Copy points to be interpolated.
        std::array<std::vector<T>,2> interpolationPoints;
        interpolationPoints[0].push_back(controlPoints[0][0]);
        interpolationPoints[1].push_back(controlPoints[1][0]);
        for (std::size_t iPoint=3ul; iPoint<controlPoints.front().size(); iPoint+=3) {
            interpolationPoints[0].push_back(controlPoints[0][iPoint]);
            interpolationPoints[1].push_back(controlPoints[1][iPoint]);
        }

        std::array<std::span<const T>,D> view;
        for (unsigned d=0u; d<D; ++d) view[d] = interpolationPoints[d];
        const auto curveDescription = geo::interpolateWithBSplineCurve<T>(view, 3);

        std::array<std::span<const T>,D> spans;
        for (unsigned i=0u; i<D; ++i)
            spans[i] = std::span<const T>(
                curveDescription[i].data(),
                curveDescription[i].size());
        std::span<const T> k(
            curveDescription.back().data(),
            curveDescription.back().size());

        return fem::maths::BSpline<T,1,D>(spans, k);
        CIE_END_EXCEPTION_TRACING
    }
}


template <class T>
void makeLineSegments(std::span<const std::vector<T>,2u> points,
                      Ref<SVG<T>> rSVG) {
    // Sanity checks.
    CIE_CHECK(
        points.size() == 2,
        std::format(
            "unexpected number of dimensions ({}) for a 2D line segment",
            points.size()))
    for (const auto& rCoordinates : points) CIE_CHECK(
        1 < rCoordinates.size(),
        std::format(
            "unexpected number of coordinates ({}) for the definition of line segment(s)",
            rCoordinates.size()))

    CIE_BEGIN_EXCEPTION_TRACING
    for (std::size_t iSegment=0ul; iSegment<points.front().size() - 1; ++iSegment) {
        std::array<std::array<T,2>,2> transformed {
            std::array<T,2> {points.front()[iSegment], points.back()[iSegment]},
            std::array<T,2> {points.front()[iSegment + 1], points.back()[iSegment + 1]}};
        rSVG.add(maths::AffineEmbedding<T,1u,2u>(transformed));
    }
    CIE_END_EXCEPTION_TRACING
}


template <class T>
void parsePathElement(Ref<const std::unordered_map<std::string,std::string>> rAttributes,
                      Ref<const maths::AffineTransform<T,2>> rTransform,
                      Ref<SVG<T>> rSVG) {
    CIE_BEGIN_EXCEPTION_TRACING
    constexpr unsigned D = 2u;

    // Find transform data.
    std::array<T,2> translation;
    {
        std::fill_n(
            translation.data(),
            D,
            static_cast<T>(0));
        const auto itTransformAttribute = rAttributes.find("transform");
        if (itTransformAttribute != rAttributes.end()) {
            const std::string_view attributeValue = itTransformAttribute->second;
            auto it = attributeValue.begin();
            const auto itEnd = attributeValue.end();
            const auto iArgumentBegin = attributeValue.find_first_of('(');
            CIE_CHECK(
                iArgumentBegin != attributeValue.npos,
                "cannot find arguments of attribute \"transform\" in element <path>")

            if (std::string_view(it, it + iArgumentBegin) == "translate") {
                it += iArgumentBegin + 1;
                for (unsigned d=0; d<D; ++d) {
                    const auto iBegin = std::string_view(it, itEnd).find_first_of("0123456789.-+Ee");
                    CIE_CHECK(
                        iBegin != std::string_view::npos,
                        std::format(
                            "cannot find argument {} of attribute \"transform\" in '{}'",
                            d, attributeValue))
                    it += iBegin;
                    auto iDelimiter = std::string_view(it, itEnd).find_first_of(", )");
                    CIE_CHECK(
                        iDelimiter != std::string_view::npos,
                        std::format(
                            "cannot find delimiter to argument {} of attribute \"transform\" in '{}'",
                            d, attributeValue))
                    // For some reason, Homebrew Clang refuses to do std::from_chars on floating point types.
                    //[[maybe_unused]] auto [pEnd, error] = std::from_chars(
                    //    it, it + iDelimiter,
                    //    translation[d]);
                    //CIE_CHECK(
                    //    error == std::errc {} && pEnd == it + iDelimiter,
                    //    std::format(
                    //        "failed to interpret component '{}' as argument {} of attribute \"transform\"",
                    //        std::string_view(it, it + iDelimiter),
                    //        d))
                    Ptr<std::string_view::value_type> pEnd = nullptr;
                    translation[d] = std::strtod(it, &pEnd);
                    CIE_CHECK(
                        pEnd == it + iDelimiter,
                        std::format(
                            "failed to interpret component '{}' as argument {} of attribute \"transform\"",
                            std::string_view(it, it + iDelimiter),
                            d))
                    it += iDelimiter;
                } // for d in range(D)
            } /*if transform is translate*/ else {
                CIE_THROW(
                    Exception,
                    std::format(
                        "unknown transform type \"{}\"",
                        std::string_view(it, it + iArgumentBegin)))
            }
        } // if itTransformAttribute != rAttributes.end()
    }

    // Expecting data defining the curve in the following format:
    // - character 'M' followed by D floats
    //      - define the offset of the curve's beginning
    // - character 'C' followed by N * 3 * D floats, defining
    //   N cubic bezier segments.
    //      - D components of the first (potentially) non-interpolated control point,
    //      - D components of the second (potentially) non-interpolated control point,
    //      - D components of the last interpolated control point.
    const auto itData = rAttributes.find("d");
    if (itData == rAttributes.end()) return;
    const std::string_view data = itData->second;
    std::array<std::vector<T>,D> points;
    char control = '\0'; // <== indicates what's being parsed in "d"
    auto it = data.begin();
    const auto itEnd = data.end();

    CIE_BEGIN_EXCEPTION_TRACING
        while (it != itEnd) {
            if (*it == 'Z') {
                // End of data. Construct whatever was being defined last.
                if (control == 'C') {
                    rSVG.add(makeBSpline<T,D>(points));
                } else if (control == 'L') {
                    makeLineSegments<T>(points, rSVG);
                }
                control = 'Z';
                break;
            } else if (*it == 'C') {
                // End of the definition of the current object.
                // Construct whatever was being defined and erase everything
                // apart from the last point.
                if (control == 'C') {
                    rSVG.add(makeBSpline<T,D>(points));
                } else if (control == 'L') {
                    makeLineSegments<T>(points, rSVG);
                }
                for (auto& rCoordinates : points) {
                    std::swap(rCoordinates.front(), rCoordinates.back());
                    rCoordinates.resize(1);
                }
                control = 'C';
                ++it;
            } else if (*it == 'M') {
                // End of the definition of the current object.
                // Construct whatever was being defined and erase everything.
                if (control == 'C') {
                    rSVG.add(makeBSpline<T,D>(points));
                } else if (control == 'L') {
                    makeLineSegments<T>(points, rSVG);
                }
                for (auto& rCoordinates : points)
                    rCoordinates.clear();
                control = 'M';
                ++it;
            } else if (*it == 'L') {
                // End of the definition of the current object.
                // Construct whatever was being defined and erase everything
                // apart from the last point.
                if (control == 'C') {
                    rSVG.add(makeBSpline<T,D>(points));
                } else if (control == 'L') {
                    makeLineSegments<T>(points, rSVG);
                }
                for (auto& rCoordinates : points) {
                    std::swap(rCoordinates.front(), rCoordinates.back());
                    rCoordinates.resize(1);
                }
                control = 'L';
                ++it;
            } else if (std::string_view(it, 1).find_first_of("0123456789.-+Ee") != std::string_view::npos) {
                // A sequence defining a pair of coordinates.
                // => Parse and store them.
                for (unsigned d=0u; d<D; ++d) {
                    // Find the next delimiter.
                    const std::size_t iComponentEnd = std::string_view(
                        it, itEnd).find_first_not_of("0123456789.-+Ee");
                    CIE_CHECK(
                        iComponentEnd != std::string_view::npos && iComponentEnd != 0ul,
                        std::format(
                            "cannot find component {} of a point in a path definition\n'{}'",
                            d, std::string_view(it, itEnd)))

                    // Parse and register the control point component.
                    const T component = std::stold(std::string(it, it + iComponentEnd));
                    points[d].push_back(component);
                    points[d].back() += translation[d];
                    it += iComponentEnd + 1;
                } // for d in range(D)
                std::array<T,D> point, transformed;
                for (unsigned d=0u; d<D; ++d) point[d] = points[d].back();
                rTransform.evaluate(point, transformed);
                for (unsigned d=0u; d<D; ++d) points[d].back() = transformed[d];
                --it;
            } else if (*it == ' ') {
                ++it;
            } else {
                CIE_THROW(
                    Exception,
                    std::format(
                        "unknown control character '{}' at index {} while parsing an SVG path\n'{}'",
                        std::distance(data.begin(), it),
                        *it,
                        data))
            }
        } // while true

        CIE_CHECK(
            control == 'Z',
            std::format(
                "failed to parse path \"d\" sequence to its end:\n{}",
                data))
    CIE_END_EXCEPTION_TRACING
    CIE_END_EXCEPTION_TRACING
}


template <class T>
void GraphML::Deserializer<fem::SVG<T>>::onElementEnd(
        Ptr<void> pThis,
        std::string_view elementName) {
    CIE_BEGIN_EXCEPTION_TRACING
    Ref<Deserializer> rThis = *static_cast<Ptr<Deserializer>>(pThis);

    if (elementName == "path") {
        parsePathElement<T>(
            rThis._attributes,
            rThis._transform,
            rThis.instance());
    } else if (elementName == "svg") {
    } else if (elementName == "defs") {
    } else {
        CIE_THROW(
            Exception,
            std::format(
                "found unknown element <{}> while parsing SVG",
                elementName))
    }

    rThis.template release<Deserializer>(&rThis, elementName);
    CIE_END_EXCEPTION_TRACING
}

template GraphML::Deserializer<fem::SVG<double>>::Deserializer(
        Ref<fem::SVG<double>>,
        Ref<GraphML::SAXHandler>,
        Ref<const maths::AffineTransform<double,2>> rTransform,
        std::string_view);

template void GraphML::Deserializer<fem::SVG<double>>::onElementBegin(
    Ptr<void> pThis,
    std::string_view elementName,
    std::span<GraphML::AttributePair> attributes);

template void GraphML::Deserializer<fem::SVG<double>>::onText(
    Ptr<void> pThis,
    std::string_view data);

template void GraphML::Deserializer<fem::SVG<double>>::onElementEnd(
    Ptr<void> pThis,
    std::string_view elementName);


} // namespace cie::fem::io
