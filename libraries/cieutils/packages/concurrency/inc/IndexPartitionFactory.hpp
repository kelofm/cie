#ifndef CIE_UTILS_INDEX_PARTITION_FACTORY_HPP
#define CIE_UTILS_INDEX_PARTITION_FACTORY_HPP

// --- Utility Includes ---
#include "packages/types/inc/types.hpp"

// --- STL Includes ---
#include <utility> // pair
#include <memory> // pair


namespace cie::mp {


struct IndexPartitionFactory
{
    using Value = Size;

    struct Partition
    {
        Value begin;

        Value end;

        long step;
    }; // struct Partition

    virtual Size size() const noexcept = 0;

    bool empty() const noexcept
    {return this->size() == 0ul;}

    virtual Partition getPartition(Size index) const = 0;

    virtual ~IndexPartitionFactory() = default;
}; // struct IndexPartitionFactory


class DynamicIndexPartitionFactory final : public IndexPartitionFactory
{
public:
    using IndexPartitionFactory::Value;

    using IndexPartitionFactory::Partition;

    DynamicIndexPartitionFactory() noexcept;

    DynamicIndexPartitionFactory(Partition range,
                                 Size partitionCount);

    ~DynamicIndexPartitionFactory() override;

    Size size() const noexcept override;

    Partition getPartition(Size index) const override;

private:
    struct Impl;
    std::unique_ptr<Impl> _pImpl;
}; // class DynamicIndexPartitionFactory


} // namespace cie::mp

#endif
