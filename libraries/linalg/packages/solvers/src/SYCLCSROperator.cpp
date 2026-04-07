// --- Linalg Includes ---
#include "packages/solvers/inc/SYCLCSROperator.hpp"

// --- Utility Includes ---
#include "hipSYCL/sycl/libkernel/group_functions.hpp"
#include "packages/macros/inc/checks.hpp"


namespace cie::linalg {


template <class TI, class TV, class TMV>
SYCLCSROperator<TI,TV,TMV>::SYCLCSROperator(
    CSRView<const TMV,const TI> lhs,
    std::shared_ptr<SYCLSpace<TV>> pSpace)
        :   _lhs(lhs),
            _pSpace(pSpace),
            _subGroupSize(0),
            _groupSize(0) {
    CIE_BEGIN_EXCEPTION_TRACING
        Ref<sycl::queue> rQueue = *_pSpace->getQueue();
        const auto subGroupSizes = rQueue.get_device().get_info<sycl::info::device::sub_group_sizes>();
        CIE_CHECK(!subGroupSizes.empty(), "")

        _subGroupSize = *std::max_element(subGroupSizes.begin(), subGroupSizes.end());
        _groupSize = (rQueue.get_device().get_info<sycl::info::device::max_work_group_size>() / _subGroupSize) * _subGroupSize;
    CIE_END_EXCEPTION_TRACING
}


template <class TI, class TV, class TMV>
void SYCLCSROperator<TI,TV,TMV>::product(
    typename Space::Value inScale,
    typename Space::ConstVectorView in,
    typename Space::Value outScale,
    typename Space::VectorView out) {
        // Sanity checks.
        CIE_CHECK(
            _lhs.rowExtents().size() - 1 == out.size() && in.size() == static_cast<std::size_t>(_lhs.columnCount()),
            std::format(
                "Incompatible matrix-vector product: [{}x{}] * [{}] = [{}]",
                _lhs.rowExtents().size() - 1, _lhs.columnCount(), in.size(), out.size()))

        CIE_BEGIN_EXCEPTION_TRACING
            Ref<sycl::queue> rQueue = *_pSpace->getQueue();
            const std::size_t itemCount = _lhs.rowCount() * _subGroupSize;
            const std::size_t itemsPerGroup = ((itemCount + _groupSize - 1) / _groupSize) * _groupSize;
            const sycl::nd_range<1> range(itemsPerGroup, _groupSize);

            rQueue.submit([&, this] (Ref<sycl::handler> rHandler) {
                Ptr<const TV> pInBegin = in.get();
                Ptr<TV> pOutBegin = out.get();
                sycl::accessor<
                    TV,
                    1,
                    sycl::access_mode::read_write,
                    sycl::access::target::local
                > groupMemory(_groupSize, rHandler);

                rHandler.parallel_for(
                    range,
                    [=, subGroupSize = _subGroupSize, lhs = _lhs] (sycl::nd_item<1> it) -> void {
                        const std::size_t iItem = it.get_global_linear_id();
                        const std::size_t iSubGroupItem = it.get_local_linear_id();
                        const std::size_t iLane = iSubGroupItem % subGroupSize;
                        const std::size_t iSubGroup = iItem / subGroupSize;
                        const std::size_t iRow = iSubGroup;

                        auto subGroup = it.get_sub_group();

                        if (iRow < lhs.rowCount()) {
                            const TI iEntryBegin = lhs.rowExtents()[iRow];
                            const TI iEntryEnd = lhs.rowExtents()[iRow + 1];

                            TV contribution = 0;
                            for (TI iEntry=iEntryBegin+iLane; iEntry<iEntryEnd; iEntry+=static_cast<TI>(subGroupSize)) {
                                const TI iColumn = lhs.columnIndices()[iEntry];
                                const TV entry = lhs.entries()[iEntry];
                                contribution += entry * pInBegin[iColumn];
                            } // for iEntry in range(iEntryBegin, iEntryEnd, subGroupSize)

                            contribution = sycl::reduce_over_group(
                                subGroup,
                                contribution,
                                sycl::plus<TV>());

                            if (subGroup.leader()) {
                                //sycl::atomic_ref<
                                //    TV,
                                //    sycl::memory_order::relaxed,
                                //    sycl::memory_scope::device,
                                //    sycl::access::address_space::global_space
                                //> ref(pOutBegin[iRow]);
                                Ref<TV> ref = pOutBegin[iRow];
                                ref *= inScale;
                                ref += outScale * contribution;
                            } // if subGroup.leader()
                        } // if iRow < lhs.rowCount()
                    }); // rHandler.parallel_for
            }).wait_and_throw(); // rQueue.submit
        CIE_END_EXCEPTION_TRACING
}


template class SYCLCSROperator<int,float,float>;
template class SYCLCSROperator<std::size_t,float,float>;
template class SYCLCSROperator<int,double,double>;
template class SYCLCSROperator<std::size_t,double,double>;


} // namespace cie::linalg
