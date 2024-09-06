#ifndef CIE_FEM_PARTITION_IMPL_HPP
#define CIE_FEM_PARTITION_IMPL_HPP

// --- FEM Includes ---
#include "packages/graph/inc/Partition.hpp"

// --- STL Includes ---
#include <limits>


namespace cie::fem {


template <class TVA, class TPA>
Partition<TVA,TPA>::Partition(RightRef<TVA> rVertexAttributes,
                              RightRef<TPA> rPolytopeAttributes) noexcept
    : Partition("",
                std::numeric_limits<Size>::max(),
                std::move(rVertexAttributes),
                std::move(rPolytopeAttributes))
{
}


template <class TVA, class TPA>
Partition<TVA,TPA>::Partition(RightRef<std::string> rName,
                              Size id,
                              RightRef<TVA> rVertexAttributes,
                              RightRef<TPA> rPolytopeAttributes) noexcept
    : PartitionBase(std::move(rName), id),
      _vertexAttributes(std::move(rVertexAttributes)),
      _polytopeAttributes(std::move(rPolytopeAttributes))
{
}


} // namespace cie::fem


#endif
