#ifndef CIE_UTILS_PARALLEL_FOR_HPP
#define CIE_UTILS_PARALLEL_FOR_HPP

// --- Internal Includes ---
#include "packages/stl_extension/inc/DynamicArray.hpp"
#include "packages/concurrency/inc/ThreadStorage.hpp"
#include "packages/concurrency/inc/ThreadPool.hpp"
#include "packages/concurrency/inc/IndexPartitionFactory.hpp"
#include "packages/compile_time/packages/concepts/inc/container_concepts.hpp"
#include "packages/types/inc/types.hpp"

// --- STL Includes ---
#include <functional>
#include <optional>


namespace cie::mp {


/** @brief Parallel for construct with thread local storage.
 *
 *  @tparam TIndex: Index type used during partitioning and the loops.
 *  @tparam TStorage: @ref ThreadStorage for storing data local to each thread.
 *  @ingroup cieutils
 */
template <concepts::Integer TIndex = Size,
          class TStorage = ThreadStorage<>>
class ParallelFor final
{
private:
    using Storage = TStorage;

    using IndexLoopFunction = typename Storage::template IndexLoopFunction<TIndex>;

    using ConstIndexLoopFunction = typename Storage::template ConstIndexLoopFunction<TIndex>;

    template <class T>
    using ObjectLoopFunction = typename Storage::template ObjectLoopFunction<T>;

    template <class T>
    using ConstObjectLoopFunction = typename Storage::template ConstObjectLoopFunction<T>;

    using Pool = ThreadPool<Storage>;

public:
    ParallelFor(ParallelFor&& rRhs) noexcept = default;

    /** @brief Construct a parallel for object with thread local storage.
     *
     *  @tparam TArgs: parameter pack of thread local storage types.
     *  @param rPool: thread pool to use for execution.
     *  @param rArgs: initializers for the thread local storages.
     *
     *  @details Each thread's local storage is initialized with the provided values,
     *           which is semantically equivalent to @c firstPrivate in OpenMP.
     *           The provided arguments are @b copied to each thread's local storage.
     *           Target functions must accept all arguments in matching order
     *           (after the index or container item).
     */
    template <class ...TArgs>
    requires (std::is_same_v<std::tuple<typename std::remove_reference_t<std::decay_t<TArgs>>...>,typename TStorage::StorageType>)
    ParallelFor(Ref<ThreadPoolBase> rPool, TArgs&&... rArgs);

    /** @brief Create a parallel for construct with initialized thread-local variables.
     *
     *  @details Each thread's local storage is initialized with the provided values,
     *           which is semantically equivalent to @c firstPrivate in OpenMP.
     *           The provided arguments are @b copied to each thread's local storage.
     *           Target functions must accept all arguments in matching order
     *           (after the index or container item).
     *  @note This is a factory function to help with reducing template clutter.
     */
    template <class ...TArgs>
    ParallelFor<TIndex,ThreadStorage<typename std::remove_reference_t<TArgs>...>>
    firstPrivate(TArgs&&... rArgs);

    /** @brief Execute an indexed for loop.
     *  @tparam TFunction Function to execute at each iteration. Must be callable
     *          with the index type followed by thread local storage types as
     *          arguments.
     *  @param rPartitionFactory factory that provides index partitions to loop over.
     *  @param rFunction target function to execute at each iteration.
     */
    template <class TFunction>
    ParallelFor& operator()(Ref<const IndexPartitionFactory> rPartitionFactory,
                            Ref<const TFunction> rFunction);

    /** @brief Execute an indexed for loop.
     *  @tparam TFunction Function to execute at each iteration. Must be callable
     *          with the index type followed by thread local storage types as
     *          arguments.
     *  @param indexMin index begin.
     *  @param indexMax index end.
     *  @param stepSize loop index increment value.
     *  @param rFunction target function to execute at each iteration.
     */
    template <class TFunction>
    ParallelFor& operator()(TIndex indexMin,
                            TIndex indexMax,
                            long stepSize,
                            Ref<const TFunction> rFunction);

    /** @brief Execute an indexed for loop beginning from 0.
     *  @tparam TFunction Function to execute at each iteration. Must be callable
     *          with the index type followed by thread local storage types as
     *          arguments.
     *  @param indexMax index end.
     *  @param rFunction target function to execute at each iteration.
     */
    template <class TFunction>
    ParallelFor& operator()(TIndex indexMax,
                            Ref<const TFunction> rFunction);

    /** @brief Execute a range-based for loop over a container.
     *  @tparam TFunction Function to execute at each iteration. Must be callable
     *          with the container's value type followed by thread local storage
     *          types as arguments.
     *  @param rContainer container holding the items to loop through.
     *  @param rFunction target function to execute at each iteration.
     */
    template <concepts::Container TContainer, class TFunction>
    ParallelFor& operator()(Ref<TContainer> rContainer,
                            Ref<const TFunction> rFunction);

    /** @brief Execute a range-based for loop over a const container.
     *  @tparam TFunction Function to execute at each iteration. Must be callable
     *          with the container's value type followed by thread local storage
     *          types as arguments.
     *  @param rContainer container holding the items to loop through.
     *  @param rFunction target function to execute at each iteration.
     */
    template <concepts::Container TContainer, class TFunction>
    ParallelFor& operator()(Ref<const TContainer> rContainer,
                            Ref<const TFunction> rFunction);

    /** @brief Execute a range-based for loop over a range.
     *  @tparam TFunction Function to execute at each iteration. Must be callable
     *          with the container's value type followed by thread local storage
     *          types as arguments.
     *  @param itBegin iterator pointing to the first item to loop through.
     *  @param itEnd iterator pointing past the last item to loop through.
     *  @param rFunction target function to execute at each iteration.
     */
    template <concepts::Iterator TIterator, class TFunction>
    ParallelFor& operator()(TIterator itBegin,
                            TIterator itEnd,
                            Ref<const TFunction> rFunction);

    Ref<const Pool> getPool() const noexcept;

    Ref<Pool> getPool() noexcept;

private:
    ParallelFor() = delete;

    ParallelFor& operator=(ParallelFor&& rRhs) = delete;

    template <class TFunction>
    void execute(Ref<const IndexPartitionFactory> rIndexPartitionFactory,
                 Ref<const TFunction> rFunction,
                 bool barrier);

    template <concepts::Container TContainer, class TFunction>
    void execute(Ref<const IndexPartitionFactory> rIndexPartitionFactory,
                 Ref<TContainer> rContainer,
                 Ref<const TFunction> rFunction);

    template <concepts::Container TContainer, class TFunction>
    void execute(Ref<const IndexPartitionFactory> rIndexPartitionFactory,
                 Ref<const TContainer> rContainer,
                 Ref<const TFunction> rFunction);

    template <concepts::Iterator TIterator, class TFunction>
    void execute(Ref<const IndexPartitionFactory> rIndexPartitionFactory,
                 TIterator itBegin,
                 Ref<const TFunction> rFunction);

private:
    Pool _pool;
};


} // namespace cie::mp

#include "packages/concurrency/impl/ParallelFor_impl.hpp"

#endif
