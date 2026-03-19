#pragma once

// --- FEM Includes ---
#include "packages/numeric/inc/CellBase.hpp"

// --- STL Includes ---
#include <tuple>
#include <span>
#include <vector>


namespace cie::fem {


template <class T>
concept DiscretizationLike
=   GraphLike<T>
&&  CellLike<typename T::Vertex::Data>
&&  maths::Expression<typename T::Data::Ansatz>
&&  maths::Expression<typename T::Data::Ansatz::Derivative>
&& requires (const typename T::Data& d) {
    {d.ansatzCount()}                       -> std::same_as<std::size_t>;
    {d.ansatz(std::size_t())}               -> std::same_as<Ref<const typename T::Data::Ansatz>>;
    {d.ansatzDerivative(std::size_t())}     -> std::same_as<Ref<const typename T::Data::Ansatz::Derivative>>;
}; // concept DiscretizationLike


template <maths::Expression TAnsatzSpace>
class MeshBase {
public:
    using Ansatz = TAnsatzSpace;

    constexpr MeshBase() noexcept = default;

    MeshBase(std::span<const TAnsatzSpace> spaces);

    [[nodiscard]] constexpr std::size_t ansatzCount() const noexcept;

    [[nodiscard]] Ref<const TAnsatzSpace> ansatz(std::size_t iAnsatz) const noexcept;

    [[nodiscard]] Ref<const typename TAnsatzSpace::Derivative> ansatzDerivative(std::size_t iAnsatz) const noexcept;

private:
    std::vector<std::tuple<
        TAnsatzSpace,
        typename TAnsatzSpace::Derivative
    >> _spaces;
}; // class MeshBase


template <maths::Expression TAnsatzSpace>
struct io::GraphML::Serializer<MeshBase<TAnsatzSpace>> {
    void header(Ref<io::GraphML::XMLElement> rElement) const;

    void operator()(
        Ref<io::GraphML::XMLElement> rElement,
        Ref<const MeshBase<TAnsatzSpace>> rInstance) const;
}; // struct Serializer


template <maths::Expression TAnsatzSpace>
struct io::GraphML::Deserializer<MeshBase<TAnsatzSpace>>
    : public io::GraphML::DeserializerBase<MeshBase<TAnsatzSpace>> {
        using io::GraphML::DeserializerBase<MeshBase<TAnsatzSpace>>::DeserializerBase;

        static void onElementBegin(
            Ptr<void> pThis,
            std::string_view elementName,
            std::span<io::GraphML::AttributePair> attributes);

        static void onText(
            Ptr<void> pThis,
            std::string_view elementName);

        static void onElementEnd(
            Ptr<const void> pThis,
            std::string_view elementName);
}; // struct Deserializer


} // namespace cie::fem

#include "packages/numeric/impl/MeshBase_impl.hpp"
