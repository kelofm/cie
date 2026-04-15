// --- Internal Includes ---
#include "embeddedPoisson2D/CellData.hpp"


namespace cie::fem {


CellData::CellData(
    VertexID id,
    AnsatzID ansatzID,
    Scalar diffusivity,
    OrientedAxes<Dimension> axes,
    RightRef<SpatialTransform> rSpatialTransform) noexcept
    : BoxBase(),
      CellBase(
        id,
        ansatzID,
        axes,
        std::move(rSpatialTransform),
        std::move(diffusivity))
{}


bool CellData::at(geo::BoxBoundable<Dimension,Scalar>::Point point) const {
    const utils::Comparison<Scalar> comparison(1e-8, 1e-10);
    StaticArray<ParametricCoordinate<Scalar>,Dimension> local;
    std::vector<Scalar> buffer(this->makeSpatialTransform().bufferSize());

    this->transform(
        Kernel<Dimension,Scalar>::cast<PhysicalCoordinate<Scalar>>(std::span<const Scalar,Dimension>(point.data(), Dimension)),
        Kernel<Dimension,Scalar>::view(local),
        buffer);

    return std::all_of(
        local.begin(),
        local.end(),
        [&comparison](Scalar coordinate) {
            return comparison.less(std::abs(coordinate), static_cast<Scalar>(1))
                || comparison.equal(std::abs(coordinate), static_cast<Scalar>(1));}
    );
}


void CellData::computeBoundingBoxImpl(Ref<BoundingBox> rBox) noexcept {
    BoundingBox::Point opposite;
    StaticArray<ParametricCoordinate<Scalar>,Dimension> localCorner;
    StaticArray<PhysicalCoordinate<Scalar>,Dimension> globalCorner;
    std::vector<Scalar> buffer(this->makeSpatialTransform().bufferSize());

    std::fill_n(
        rBox.base().data(),
        Dimension,
        std::numeric_limits<Scalar>::max());
    std::fill_n(
        opposite.data(),
        Dimension,
        std::numeric_limits<Scalar>::lowest());

    StaticArray<std::uint8_t,Dimension> state {0u, 0u};
    StaticArray<ParametricCoordinate<Scalar>,2> ordinates {-1.0, 1.0};

    do {
        // Compute the corner in local space.
        std::transform(
        state.begin(),
            state.end(),
            localCorner.begin(),
            [&ordinates](std::uint8_t iOrdinate) {
                return ordinates[iOrdinate];
            });

        // Transform the corner to global space.
        this->transform(
            Kernel<Dimension,Scalar>::view(localCorner),
            Kernel<Dimension,Scalar>::view(globalCorner),
            buffer);

        // Extend box definition.
        for (unsigned iDimension=0u; iDimension<Dimension; ++iDimension) {
            rBox.base()[iDimension] = std::min<Scalar>(rBox.base()[iDimension], globalCorner[iDimension]);
            opposite[iDimension] = std::max<Scalar>(opposite[iDimension], globalCorner[iDimension]);
        } // for iDimension in range(Dimension)
    } while (cie::maths::OuterProduct<Dimension>::next(2u, state.data()));

    std::transform(
        opposite.begin(),
        opposite.end(),
        rBox.base().begin(),
        rBox.lengths().begin(),
        std::minus<Scalar>());
}


} // namespace cie::fem
