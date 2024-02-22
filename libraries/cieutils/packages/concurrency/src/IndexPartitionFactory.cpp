// --- Utility Includes ---
#include "packages/concurrency/inc/IndexPartitionFactory.hpp"
#include "packages/macros/inc/checks.hpp"
#include "packages/stl_extension/inc/DynamicArray.hpp"


namespace cie::mp {


struct DynamicIndexPartitionFactory::Impl
{
    DynamicArray<Size> partitions;

    long step;
}; // DynamicIndexPartitionFactory::Impl


DynamicIndexPartitionFactory::DynamicIndexPartitionFactory() noexcept
    : _pImpl(new Impl)
{
}


DynamicIndexPartitionFactory::DynamicIndexPartitionFactory(Partition range,
                                                           Size partitionCount)
    : DynamicIndexPartitionFactory()
{
    CIE_BEGIN_EXCEPTION_TRACING

    // Check index range
    CIE_CHECK(partitionCount != 0, "cannot generate index partition for 0 threads")
    CIE_CHECK(range.step != 0, "cannot generate index partition with a step size of 0")
    CIE_CHECK((range.begin <= range.end) == (0 < range.step), "step direction is inconsistent with range definition")

    // Handle reverse loops
    if (range.step < 0) {
        std::swap(range.begin, range.end);
        range.step = -range.step;
    }

    // Initialize index blocks
    Value blockSize = (range.end - range.begin) / range.step / partitionCount;
    blockSize += blockSize % range.step;
    blockSize = std::max(blockSize, 1ul);

    _pImpl->partitions.reserve(partitionCount + 1);
    Value tmp = range.begin;
    for (Size i=0; i<partitionCount && tmp<range.end; ++i, tmp+=blockSize) {
        _pImpl->partitions.push_back(tmp);
    }
    _pImpl->partitions.push_back(range.end);
    _pImpl->step = range.step;

    CIE_END_EXCEPTION_TRACING
}


DynamicIndexPartitionFactory::~DynamicIndexPartitionFactory()
{
}


Size DynamicIndexPartitionFactory::size() const noexcept
{
    return _pImpl->partitions.empty() ? 0ul : _pImpl->partitions.size() - 1;
}


DynamicIndexPartitionFactory::Partition DynamicIndexPartitionFactory::getPartition(Size index) const
{
    CIE_OUT_OF_RANGE_CHECK(std::max(index, index + 1) < _pImpl->partitions.size())
    return {_pImpl->partitions[index], _pImpl->partitions[index + 1], _pImpl->step};
}


} // namespace cie::mp
