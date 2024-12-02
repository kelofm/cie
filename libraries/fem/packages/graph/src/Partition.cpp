// --- FEM Includes ---
#include "packages/graph/inc/Partition.hpp"

// --- STL Includes ---
#include <limits>


namespace cie::fem {


PartitionBase::PartitionBase() noexcept
    : PartitionBase("", std::numeric_limits<Size>::max())
{
}


PartitionBase::PartitionBase(RightRef<std::string> rName, PartitionID id)
    : utils::NamedObject(std::move(rName)),
      utils::IDObject<PartitionID>(id)
{
}


} // namespace cie::fem
