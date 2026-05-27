#ifdef CIE_ENABLE_SYCL

// --- Linalg Includes ---
#include "packages/solvers/inc/SYCLCSRSubsetOperator.hpp"

// --- Utility Includes ---
#include "packages/macros/inc/checks.hpp"


namespace cie::linalg {


template <class TI, class TV>
SYCLCSRSubsetOperator<TI,TV>::SYCLCSRSubsetOperator(
    CSRView<const TV,const TI> lhs,
    typename IndexSpace::ConstVectorView subset,
    std::shared_ptr<Space> pSpace,
    std::shared_ptr<IndexSpace> pIndexSpace)
        :   _lhs(lhs),
            _subset(subset),
            _pSpace(pSpace),
            _pIndexSpace(pIndexSpace),
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


template <class TI, class TV>
void SYCLCSRSubsetOperator<TI,TV>::product(
    typename Space::Value inScale,
    typename Space::ConstVectorView in,
    typename Space::Value outScale,
    typename Space::VectorView out) {
        // Sanity checks.
        CIE_CHECK(
            _lhs.rowCount() == out.size() && in.size() == static_cast<std::size_t>(_lhs.columnCount()),
            std::format(
                "Incompatible matrix-vector product: [{}x{}] * [{}] = [{}]",
                _lhs.rowCount(), _lhs.columnCount(), in.size(), out.size()))

        CIE_BEGIN_EXCEPTION_TRACING
            Ref<sycl::queue> rQueue = *_pSpace->getQueue();
            const std::size_t itemCount = _lhs.rowCount() * _subGroupSize;
            const std::size_t itemsPerGroup = ((itemCount + _groupSize - 1) / _groupSize) * _groupSize;
            const sycl::nd_range<1> range(itemsPerGroup, _groupSize);

            rQueue.submit([&, this] (Ref<sycl::handler> rHandler) {
                Ptr<const TV> pInBegin = in.get();
                Ptr<TV> pOutBegin = out.get();

                rHandler.parallel_for(
                    range,
                    [inScale, pInBegin, outScale, pOutBegin, subGroupSize = _subGroupSize, lhs = _lhs] (sycl::nd_item<1> it) -> void {
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
                                pOutBegin[iRow] = inScale * pOutBegin[iRow] + outScale * contribution;
                            } // if subGroup.leader()
                        } // if iRow < lhs.rowCount()
                    }); // rHandler.parallel_for
            }).wait_and_throw(); // rQueue.submit
        CIE_END_EXCEPTION_TRACING
}


template class SYCLCSRSubsetOperator<int,float>;
template class SYCLCSRSubsetOperator<std::size_t,float>;
template class SYCLCSRSubsetOperator<int,double>;
template class SYCLCSRSubsetOperator<std::size_t,double>;


} // namespace cie::linalg

#endif
