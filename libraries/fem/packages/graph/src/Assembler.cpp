// --- FEM Includes ---
#include "packages/graph/inc/Assembler.hpp"

// --- Utility Includes ---
#include "packages/concurrency/inc/Mutex.hpp"
#include "packages/concurrency/inc/ParallelFor.hpp"

// --- STL Includes ---
#include <numeric> // inclusive_scan


namespace cie::fem {


Assembler::Assembler() noexcept
    : _dofCounter(0)
{}


Assembler::Assembler(std::size_t dofBegin) noexcept
    : _dofCounter(dofBegin)
{}


std::size_t Assembler::dofCount() const noexcept {
    return _dofCounter;
}


template <concepts::Integer T>
void Assembler::reorder(std::span<const T> map) {
    CIE_CHECK(_dofCounter <= map.size(), "")
    CIE_BEGIN_EXCEPTION_TRACING
        for (auto it=_dofMap.begin(); it!=_dofMap.end(); ++it) {
            for (auto& rDof : it.value())
                if (rDof.has_value()) {
                    rDof = static_cast<std::size_t>(map[*rDof]);
                }
        }
    CIE_END_EXCEPTION_TRACING
}


#define CIE_INSTANTIATE_REORDER(T)                          \
    template void Assembler::reorder<T>(std::span<const T>);


CIE_INSTANTIATE_REORDER(std::size_t)
CIE_INSTANTIATE_REORDER(int)


#undef CIE_INSTANTIATE_REORDER


template <class TIndex, class TValue>
void Assembler::makeCSRMatrix(
    Ref<TIndex> rRowCount,
    Ref<TIndex> rColumnCount,
    Ref<DynamicArray<TIndex>> rRowExtents,
    Ref<DynamicArray<TIndex>> rColumnIndices,
    Ref<DynamicArray<TValue>> rNonzeros,
    OptionalRef<mp::ThreadPoolBase> rThreadPool) const {
        CIE_BEGIN_EXCEPTION_TRACING

        // Initialize matrix containers
        rRowExtents.clear();
        rColumnIndices.clear();
        rNonzeros.clear();

        const TIndex rowCount = this->dofCount();
        rRowCount = rowCount;
        rColumnCount = rowCount;

        // Initialize mutexes
        // Each mutex is assigned a range of row indices it provides access to
        const TIndex threadCount = rThreadPool.has_value() ? rThreadPool.value().size() : 1;
        const TIndex minRowsPerMutex = rowCount / threadCount;
        const TIndex leftoverRowCount = rowCount % threadCount;
        DynamicArray<TIndex> mutexExtents;
        DynamicArray<mp::Mutex<tags::SMP>> mutexes(threadCount);

        mutexExtents.reserve(threadCount + 1);
        for (TIndex iThread=0; iThread<threadCount; ++iThread) {
            mutexExtents.push_back(
                mutexExtents.empty()
                ? 0
                : mutexExtents.back()
                + minRowsPerMutex
                + (iThread < leftoverRowCount ? 1 : 0));
        }
        mutexExtents.push_back(rowCount);

        // Insert column indices into temporary row container
        DynamicArray<DynamicArray<TIndex>> columnIndices(rowCount);

        {
            const auto job = [&columnIndices, &mutexExtents, &mutexes](const auto& rIndices) -> void {
                for (const auto iRow : rIndices) {
                    // Find which thread the row belongs to.
                    const auto itMutexExtent = std::max(
                        mutexExtents.begin(),
                        std::lower_bound(
                            mutexExtents.begin(),
                            mutexExtents.end(),
                            iRow) - 1
                    );
                    CIE_OUT_OF_RANGE_CHECK(itMutexExtent < mutexExtents.end());
                    auto& rMutex = mutexes[std::distance(mutexExtents.begin(), itMutexExtent)];

                    std::scoped_lock<mp::Mutex<tags::SMP>> lock(rMutex);
                    CIE_OUT_OF_RANGE_CHECK(iRow < columnIndices.size());
                    columnIndices[iRow].insert(
                        columnIndices[iRow].end(),
                        rIndices.begin(),
                        rIndices.end());
                } // for iRow in rIndices
            };

            if (threadCount < 2) {
                for (const auto& rDofIndices : this->values()) job(rDofIndices);
            } else {
                const auto& rDofIndexContainers = this->values();
                mp::ParallelFor<>(rThreadPool.value())(
                    rDofIndexContainers.begin(),
                    rDofIndexContainers.end(),
                    job);
            }
        }

        // Make column indices sorted and unique
        {
            const auto job = [] (Ref<DynamicArray<TIndex>> rIndices) {
                std::sort(rIndices.begin(), rIndices.end());
                rIndices.erase(std::unique(rIndices.begin(), rIndices.end()), rIndices.end());
            };

            if (threadCount < 2) {
                for (auto& rIndices : columnIndices) job(rIndices);
            } else {
                mp::ParallelFor<>(rThreadPool.value())(columnIndices, job);
            }
        }

        // Compute row extents
        rRowExtents.resize(rowCount + 1);
        rRowExtents.front() = 0;

        std::inclusive_scan(
            columnIndices.begin(),
            columnIndices.end(),
            rRowExtents.begin() + 1,
            [](TIndex left, const auto& rRight) -> TIndex {return left + rRight.size();},
            static_cast<TIndex>(0));

        // Copy column indices to CSR container
        rColumnIndices.resize(rRowExtents.back());
        rNonzeros.resize(rRowExtents.back());

        {
            const auto job = [&columnIndices, &rColumnIndices, &rRowExtents] (const TIndex iRow) {
                std::copy(
                    columnIndices[iRow].begin(),
                    columnIndices[iRow].end(),
                    rColumnIndices.begin() + rRowExtents[iRow]);
            };

            if (threadCount < 2) {
                for (TIndex iRow=0; iRow<rowCount; ++iRow) job(iRow);
            } else {
                mp::ParallelFor<>(rThreadPool.value())(rowCount, job);
            }
        }

        CIE_END_EXCEPTION_TRACING
}


#define CIE_INSTANTIATE_CSR_FACTORY(TIndex, TValue)         \
    template void Assembler::makeCSRMatrix<TIndex,TValue>(  \
        Ref<TIndex>,                                        \
        Ref<TIndex>,                                        \
        Ref<DynamicArray<TIndex>>,                          \
        Ref<DynamicArray<TIndex>>,                          \
        Ref<DynamicArray<TValue>>,                          \
        OptionalRef<mp::ThreadPoolBase>) const;

CIE_INSTANTIATE_CSR_FACTORY(int, float)
CIE_INSTANTIATE_CSR_FACTORY(int, double)
CIE_INSTANTIATE_CSR_FACTORY(std::size_t, float)
CIE_INSTANTIATE_CSR_FACTORY(std::size_t, double)

#undef CIE_INSTANTIATE_CSR_FACTORY


template <unsigned D, class T>
void makeAnsatzMask(
    Ref<const Assembler> rAssembler,
    std::size_t setSize,
    std::span<T> mask) {
        CIE_CHECK(rAssembler.dofCount() <= mask.size(), "")
        std::vector<T> localMask(intPow(setSize, D));
        makeAnsatzMask<D,T>(setSize, localMask);
        for (const auto& dofs : rAssembler.values()) {
            assert(dofs.size() == mask.size());
            for (std::size_t iDoF=0ul; iDoF<dofs.size(); ++iDoF)
                mask[dofs[iDoF]] = localMask[iDoF];
        }
}


#define CIE_INSTANTIATE_MASK_FACTORY(TIndex)                                                    \
    template void makeAnsatzMask<1,TIndex>(Ref<const Assembler>,std::size_t,std::span<TIndex>); \
    template void makeAnsatzMask<2,TIndex>(Ref<const Assembler>,std::size_t,std::span<TIndex>); \
    template void makeAnsatzMask<3,TIndex>(Ref<const Assembler>,std::size_t,std::span<TIndex>);

CIE_INSTANTIATE_MASK_FACTORY(std::uint8_t)
CIE_INSTANTIATE_MASK_FACTORY(std::uint16_t)
CIE_INSTANTIATE_MASK_FACTORY(std::size_t)
CIE_INSTANTIATE_MASK_FACTORY(int)
CIE_INSTANTIATE_MASK_FACTORY(float)
CIE_INSTANTIATE_MASK_FACTORY(double)

#undef CIE_INSTANTIATE_MASK_FACTORY


} // namespace cie::fem
