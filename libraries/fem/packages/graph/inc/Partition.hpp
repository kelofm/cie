#ifndef CIE_FEM_PARTITION_HPP
#define CIE_FEM_PARTITION_HPP

// --- FEM Includes ---
#include "packages/utilities/inc/AttributeContainer.hpp"

// --- Utility Includes ---
#include "packages/stl_extension/inc/DynamicArray.hpp"
#include "packages/types/inc/NamedObject.hpp"
#include "packages/types/inc/IDObject.hpp"
#include <packages/stl_extension/inc/StrongTypeDef.hpp>

// --- STL Includes ---
#include <concepts>
#include <memory>


namespace cie::fem {
class PartitionBase;
} // namespace cie::fem


namespace cie::concepts {


template <class T>
concept Partition
= std::derived_from<T,fem::PartitionBase> && requires () {
    typename T::VertexAttributes;
    typename T::PolytopeAttributes;
};


} // namespace cie::concepts


namespace cie::fem {


class PartitionBase : public utils::NamedObject,
                      public utils::IDObject<Size>
{
public:
    enum Item
    {
        Vertex,
        Polytope
    }; // enum Item

public:
    PartitionBase() noexcept;

    PartitionBase(const PartitionBase&) = delete;

    PartitionBase(PartitionBase&&) noexcept = default;

    PartitionBase(RightRef<std::string> r_name, Size id);

    PartitionBase& operator=(const PartitionBase&) = delete;

    PartitionBase& operator=(PartitionBase&&) noexcept = default;

    virtual ~PartitionBase() = default;
}; // class PartitionBase



template <class TVertexAttributes,
          class TPolytopeAttributes>
class Partition : public PartitionBase
{
public:
    using VertexAttributes = TVertexAttributes;

    using PolytopeAttributes = TPolytopeAttributes;

    static_assert(concepts::AttributeContainer<VertexAttributes>);

    static_assert(concepts::AttributeContainer<PolytopeAttributes>);

public:
    Partition() noexcept = default;

    Partition(RightRef<VertexAttributes> r_vertexAttributes,
              RightRef<PolytopeAttributes> r_polytopeAttributes) noexcept;

    Partition(RightRef<std::string> r_name,
              Size id,
              RightRef<VertexAttributes> r_vertexAttributes,
              RightRef<PolytopeAttributes> r_polytopeAttributes) noexcept;

private:
    friend class PartitionManager;

    VertexAttributes _vertexAttributes;

    PolytopeAttributes _polytopeAttributes;
}; // class Partition



CIE_STRONG_TYPEDEF(Size, ParentIndex);



} // namespace cie::fem

#include "packages/graph/impl/Partition_impl.hpp"

#endif
