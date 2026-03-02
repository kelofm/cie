#pragma once

// --- FEM Includes ---
#include "packages/maths/inc/BSpline.hpp"
#include "packages/maths/inc/AffineEmbedding.hpp"
#include "packages/maths/inc/AffineTransform.hpp"
#include "packages/io/inc/GraphML.hpp"

// --- STL Includes ---
#include <vector>
#include <variant>
#include <unordered_map>


namespace cie::fem {


template <concepts::Numeric TValue>
class SVG {
public:
    using Value = TValue;

    using CurveVariant = std::variant<
        std::monostate,
        maths::BSpline<TValue,1u,2u>,
        maths::AffineEmbedding<TValue,1u,2u>>;

    SVG() noexcept = default;

    SVG(RightRef<std::vector<CurveVariant>> rCurves);

    void add(RightRef<CurveVariant> rCurve);

    [[nodiscard]] std::vector<std::array<TValue,2>> tesselate(std::size_t segmentsPerCurve) const;

    constexpr std::size_t curveCount() const noexcept {
        return _curves.size();
    }

private:
    std::vector<CurveVariant> _curves;
}; // class SVG


} // namespace cie::fem


namespace cie::fem::io {


template <class T>
struct GraphML::Deserializer<fem::SVG<T>>
    : public GraphML::DeserializerBase<fem::SVG<T>> {
    Deserializer(
        Ref<fem::SVG<T>> rInstance,
        Ref<SAXHandler> rSAX,
        Ref<const maths::AffineTransform<T,2>> rTransform = {},
        std::string_view parseMode = "svg");

    static void onElementBegin(
        Ptr<void> pThis,
        std::string_view elementName,
        std::span<GraphML::AttributePair> attributes);

    static void onText(
        Ptr<void> pThis,
        std::string_view data);

    static void onElementEnd(
        Ptr<void> pThis,
        std::string_view elementName);

private:
    std::string _mode;

    std::unordered_map<std::string,std::string> _attributes;

    maths::AffineTransform<T,2> _transform;
}; // struct GraphML::Deserializer<SVG>


} // namespace cie::fem::io
