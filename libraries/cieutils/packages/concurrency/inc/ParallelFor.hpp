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
    ParallelFor(ParallelFor&& r_rhs) = default;

    ParallelFor& operator=(ParallelFor&& r_rhs) = default;

    /** @brief Construct a parallel for object with thread local storage.
     *
     *  @tparam TArgs: parameter pack of thread local storage types.
     *  @param r_args: initializers for the thread local storages.
     *
     *  @details Each thread's local storage is initialized with the provided values,
     *           which is semantically equivalent to @c firstPrivate in OpenMP.
     *           The provided arguments are @b copied to each thread's local storage.
     *           Target functions must accept all arguments in matching order
     *           (after the index or container item).
     */
    template <class ...TArgs>
    requires (!std::is_same_v<std::tuple<typename std::decay<TArgs>::type...>,std::tuple<ParallelFor<TIndex,TStorage>>>)
    ParallelFor(TArgs&&... r_args);

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
    static ParallelFor<TIndex,ThreadStorage<typename std::remove_reference<TArgs>::type...>>
    firstPrivate(TArgs&&... r_args);

    /// @brief Set the ThreadPool that will be used during the for loop.
    void setPool(RightRef<Pool> r_pool) noexcept;

    /// @brief Set the ThreadPool that will be used during the for loop.
    void setPool(Ref<Pool> r_pool);

    /** @brief Execute an indexed for loop.
     *
     *  @param indexMin index begin
     *  @param indexMax index end
     *  @param stepSize loop index increment value
     *  @param r_function target function to execute at each iteration
     */
    ParallelFor& operator()(TIndex indexMin,
                            TIndex indexMax,
                            TIndex stepSize,
                            IndexLoopFunction&& r_function);

    /** @brief Execute an indexed for loop beginning from 0.
     *
     *  @param indexMax index end
     *  @param r_function target function to execute at each iteration
     */
    ParallelFor& operator()(TIndex indexMax,
                            IndexLoopFunction&& r_function);

    /** @brief Execute a range-based for loop over a container.
     *
     *  @param r_container container holding the items to loop through
     *  @param r_function target function to execute at each iteration
     */
    template <concepts::Container TContainer>
    ParallelFor& operator()(TContainer& r_container,
                            ObjectLoopFunction<typename TContainer::value_type>&& r_function);

    /** @brief Execute a range-based for loop over a const container.
     *
     *  @param r_container container holding the items to loop through
     *  @param r_function target function to execute at each iteration
     */
    template <concepts::Container TContainer>
    ParallelFor& operator()(const TContainer& r_container,
                            ConstObjectLoopFunction<typename TContainer::value_type>&& r_function);

    /** @brief Execute a range-based for loop over a container.
     *
     *  @param r_container container holding the items to loop through
     *  @param r_function target function to execute at each iteration
     */
    template <concepts::Container TContainer>
    ParallelFor& operator()(TContainer& r_container,
                            ObjectLoopFunction<typename TContainer::value_type>& r_function);

    /** @brief Execute a range-based for loop over a const container.
     *
     *  @param r_container container holding the items to loop through
     *  @param r_function target function to execute at each iteration
     */
    template <concepts::Container TContainer>
    ParallelFor& operator()(const TContainer& r_container,
                            ConstObjectLoopFunction<typename TContainer::value_type>& r_function);

    Ref<const Pool> getPool() const noexcept;

    Ref<Pool> getPool() noexcept;

protected:
    IndexPartition makeIndexPartition(TIndex indexMin,
                                      TIndex indexMax,
                                      TIndex stepSize);

    void execute(const IndexPartition& r_indexPartition,
                 TIndex stepSize,
                 IndexLoopFunction&& r_function);

    template <concepts::Container TContainer>
    void execute(const IndexPartition& r_indexPartition,
                 TIndex stepSize,
                 TContainer& r_container,
                 ObjectLoopFunction<typename TContainer::value_type>&& r_function);

private:
    Pool _pool;
};


} // namespace cie::mp

#include "packages/concurrency/impl/ParallelFor_impl.hpp"

#endif