#ifndef CIE_FEM_GRAPH_BOUNDARY_ID_IMPL_HPP
#define CIE_FEM_GRAPH_BOUNDARY_ID_IMPL_HPP

// --- Utility Includes ---
#include "packages/maths/inc/bit.hpp"

// --- FEM Includes ---
#include "packages/graph/inc/BoundaryID.hpp"

// --- STL Includes ---
#include <climits>


namespace cie::fem {


inline constexpr BoundaryID::BoundaryID() noexcept
    : BoundaryID(0, false)
{
}


inline constexpr BoundaryID::BoundaryID(const char name[3])
{
    bool direction = false;
    unsigned dimension = 0;

    switch (name[0]) {
        case '-': {direction = false; break;}
        case '+': {direction = true; break;}
        default: throw std::runtime_error("Invalid boundary name " + std::string(name));
    } // switch name[0]

    switch (name[1]) {
        case 'x': {dimension = 0u; break;}
        case 'y': {dimension = 1u; break;}
        case 'z': {dimension = 2u; break;}
        default: {
            if ('0' <= name[1] && name[1] <= '9') {
                dimension = static_cast<unsigned>(name[1] - '0');
            } else {
                throw std::runtime_error("Invalid boundary name " + std::string(name));
            }
        }
    } // switch name[1]

    *this = BoundaryID(dimension, direction);
}


inline constexpr BoundaryID::BoundaryID(unsigned dimension, bool direction)
    : _id((1u << ++dimension) + direction)
{
    #ifdef CIE_ENABLE_OUT_OF_RANGE_CHECKS
    if (CHAR_BIT * sizeof(unsigned) <= dimension || dimension == 0) {
        throw std::invalid_argument("Input dimension out of range");
    }
    #endif
}


inline constexpr BoundaryID& BoundaryID::operator++() noexcept
{
    _id <<= (1u ^ ((_id ^= 1u) & 1u));
    return *this;
}


inline constexpr BoundaryID BoundaryID::operator++(int) noexcept
{
    BoundaryID copy(*this);
    ++*this;
    return copy;
}


inline constexpr BoundaryID BoundaryID::operator-() const noexcept
{
    BoundaryID copy = *this;
    copy._id ^= 1u;
    return copy;
}


inline unsigned BoundaryID::getDimension() const noexcept
{
    return utils::getNumberOfTrailingZeros(_id & (~1u)) - 1u;
}


inline constexpr bool BoundaryID::getDirection() const noexcept
{
    return _id & 1u;
}


inline constexpr unsigned BoundaryID::getData() const noexcept
{
    return _id;
}


inline constexpr bool operator==(BoundaryID left, BoundaryID right) noexcept
{
    return left._id == right._id;
}


inline constexpr bool operator!=(BoundaryID left, BoundaryID right) noexcept
{
    return left._id != right._id;
}


inline constexpr bool operator<(BoundaryID left, BoundaryID right) noexcept
{
    return left._id < right._id;
}


inline void io::GraphML::Serializer<BoundaryID>::header(Ref<XMLElement> rElement)
{
    XMLElement defaultData = rElement.addChild("default");
    std::stringstream stream;
    stream << BoundaryID();
    defaultData.setValue(stream.view());
}


inline void io::GraphML::Serializer<BoundaryID>::operator()(Ref<XMLElement> rElement,
                                                            Ref<const BoundaryID> rInstance)
{
    std::stringstream stream;
    stream << rInstance;
    rElement.setValue(stream.view());
}


inline void io::GraphML::Deserializer<BoundaryID>::onElementBegin(Ptr<void> ,
                                                                  std::string_view,
                                                                  std::span<GraphML::AttributePair>) noexcept
{
}


inline void io::GraphML::Deserializer<BoundaryID>::onText(Ptr<void> pThis,
                                                          std::string_view data)
{
    if (data.size() != 2) {
        CIE_THROW(
            Exception,
            "Expecting a 2-character string for a BoundaryID, but got \"" << data << "\"."
        )
    }

    Ref<Deserializer> rThis = *static_cast<Ptr<Deserializer>>(pThis);
    rThis.instance() = BoundaryID(data.data());
}


inline void io::GraphML::Deserializer<BoundaryID>::onElementEnd(Ptr<void> pThis,
                                                                std::string_view elementName) noexcept
{
    Ref<Deserializer> rThis = *static_cast<Ptr<Deserializer>>(pThis);
    rThis.template release<Deserializer>(&rThis, elementName);
}


} // namespace cie::fem


namespace std {


inline size_t hash<cie::fem::BoundaryID>::operator()(cie::fem::BoundaryID id) const noexcept
{
    return hash<unsigned>()(id.getData());
}


} // namespace std


#endif
