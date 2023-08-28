// --- FEM Includes ---
#include "packages/graph/inc/Partition.hpp"

// --- STL Includes ---
#include <limits>


namespace cie::fem {


PartitionBase::PartitionBase() noexcept
    : PartitionBase("", std::numeric_limits<Size>::max())
{
}


PartitionBase::PartitionBase(RightRef<std::string> r_name, Size id)
    : utils::NamedObject(std::move(r_name)),
      utils::IDObject<Size>(id)
{
}


} // namespace cie::fem
