#ifndef CIE_UTILS_PARALLEL_FOR_HPP
#define CIE_UTILS_PARALLEL_FOR_HPP

// --- Internal Includes ---
#include "packages/stl_extension/inc/DynamicArray.hpp"
#include "packages/concurrency/inc/ThreadStorage.hpp"
#include "packages/concurrency/inc/ThreadPool.hpp"
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
    using IndexPartition = DynamicArray<TIndex>;

    using Storage = TStorage;

    using IndexLoopFunction = typename Storage::template IndexLoopFunction<TIndex>;

    using ConstIndexLoopFunction = typename Storage::template ConstIndexLoopFunction<TIndex>;

    template <class T>
    using ObjectLoopFunction = typename Storage::template ObjectLoopFunction<T>;

    template <class T>
    using ConstObjectLoopFunction = typename Storage::template ConstObjectLoopFunction<T>;

    using Pool = ThreadPool<Storage>;

public:
    ParallelFor(ParallelFor&& r_rhs) noexcept = default;

    /** @brief Construct a parallel for object with thread local storage.
     *
     *  @tparam TArgs: parameter pack of thread local storage types.
     *  @param r_pool: thread pool to use for execution.
     *  @param r_args: initializers for the thread local storages.
     *
     *  @details Each thread's local storage is initialized with the provided values,
     *           which is semantically equivalent to @c firstPrivate in OpenMP.
     *           The provided arguments are @b copied to each thread's local storage.
     *           Target functions must accept all arguments in matching order
     *           (after the index or container item).
     */
    template <class ...TArgs>
    requires (std::is_same_v<std::tuple<typename std::remove_reference_t<std::decay_t<TArgs>>...>,typename TStorage::StorageType>)
    ParallelFor(Ref<ThreadPoolBase> r_pool, TArgs&&... r_args);

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
    firstPrivate(TArgs&&... r_args);

    /** @brief Execute an indexed for loop.
     *  @tparam TFunction Function to execute at each iteration. Must be callable
     *          with the index type followed by thread local storage types as
     *          arguments.
     *  @param indexMin index begin.
     *  @param indexMax index end.
     *  @param stepSize loop index increment value.
     *  @param r_function target function to execute at each iteration.
     */
    template <class TFunction>
    ParallelFor& operator()(TIndex indexMin,
                            TIndex indexMax,
                            TIndex stepSize,
                            Ref<const TFunction> r_function);

    /** @brief Execute an indexed for loop beginning from 0.
     *  @tparam TFunction Function to execute at each iteration. Must be callable
     *          with the index type followed by thread local storage types as
     *          arguments.
     *  @param indexMax index end
     *  @param r_function target function to execute at each iteration
     */
    template <class TFunction>
    ParallelFor& operator()(TIndex indexMax,
                            Ref<const TFunction> r_function);

    /** @brief Execute a range-based for loop over a container.
     *  @tparam TFunction Function to execute at each iteration. Must be callable
     *          with the container's value type followed by thread local storage
     *          types as arguments.
     *  @param r_container container holding the items to loop through
     *  @param r_function target function to execute at each iteration
     */
    template <concepts::Container TContainer, class TFunction>
    ParallelFor& operator()(Ref<TContainer> r_container,
                            Ref<const TFunction> r_function);

    /** @brief Execute a range-based for loop over a const container.
     *  @tparam TFunction Function to execute at each iteration. Must be callable
     *          with the container's value type followed by thread local storage
     *          types as arguments.
     *  @param r_container container holding the items to loop through
     *  @param r_function target function to execute at each iteration
     */
    template <concepts::Container TContainer, class TFunction>
    ParallelFor& operator()(Ref<const TContainer> r_container,
                            Ref<const TFunction> r_function);

    Ref<const Pool> getPool() const noexcept;

    Ref<Pool> getPool() noexcept;

private:
    ParallelFor() = delete;

    ParallelFor& operator=(ParallelFor&& r_rhs) = delete;

    IndexPartition makeIndexPartition(TIndex indexMin,
                                      TIndex indexMax,
                                      TIndex stepSize);

    template <class TFunction>
    void execute(Ref<const IndexPartition> r_indexPartition,
                 TIndex stepSize,
                 Ref<const TFunction> r_function);

    template <concepts::Container TContainer, class TFunction>
    void execute(Ref<const IndexPartition> r_indexPartition,
                 TIndex stepSize,
                 Ref<TContainer> r_container,
                 Ref<const TFunction> r_function);

    template <concepts::Container TContainer, class TFunction>
    void execute(Ref<const IndexPartition> r_indexPartition,
                 TIndex stepSize,
                 Ref<const TContainer> r_container,
                 Ref<const TFunction> r_function);

private:
    Pool _pool;
};


} // namespace cie::mp

#include "packages/concurrency/impl/ParallelFor_impl.hpp"

#endif