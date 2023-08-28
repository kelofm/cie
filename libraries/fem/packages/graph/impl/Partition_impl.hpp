#ifndef CIE_FEM_PARTITION_IMPL_HPP
#define CIE_FEM_PARTITION_IMPL_HPP

// --- FEM Includes ---
#include "packages/graph/inc/Partition.hpp"

// --- STL Includes ---
#include <limits>


namespace cie::fem {


template <class TVA, class TPA>
Partition<TVA,TPA>::Partition(RightRef<TVA> r_vertexAttributes,
                              RightRef<TPA> r_polytopeAttributes) noexcept
    : Partition("",
                std::numeric_limits<Size>::max(),
                std::move(r_vertexAttributes),
                std::move(r_polytopeAttributes))
{
}


template <class TVA, class TPA>
Partition<TVA,TPA>::Partition(RightRef<std::string> r_name,
                              Size id,
                              RightRef<TVA> r_vertexAttributes,
                              RightRef<TPA> r_polytopeAttributes) noexcept
    : PartitionBase(std::move(r_name), id),
      _vertexAttributes(std::move(r_vertexAttributes)),
      _polytopeAttributes(std::move(r_polytopeAttributes))
{
}


} // namespace cie::fem


#endif
