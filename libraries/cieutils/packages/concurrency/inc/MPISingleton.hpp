#ifndef CIE_MPI_SINGLETON_HPP
#define CIE_MPI_SINGLETON_HPP

// --- Utility Includes ---
#include "packages/concurrency/inc/MPI.hpp"

// --- STL Includes ---
#include <mutex>
#include <memory>


namespace cie::mpi {


class MPISingleton
{
private:
    using Mutex = std::mutex;

    using Lock = std::scoped_lock<Mutex>;

public:
    /// @brief Check whether the static MPI interface is constructed and initialized.
    static bool isInitialized();

    /// @brief Check whether the static MPI interface is constructed and distributed.
    static bool isDistributed();

    /** @brief Get the static MPI interface.
     *
     *  @details If the static MPI interface was not constructed and initialized yet,
     *           it will be after this function is called. If CiE was not compiled
     *           wiht MPI support, a serial dummy MPI interface is constructed.
     *  @note If you need to control what kind of MPI interface to construct, use
     *        @ref MPISingleton::get(Ptr<MPI_Comm>).
     */
    static Ref<MPI> get();

    /** @brief Get the static MPI interface.
     *
     *  @param p_comm: Ptr to the MPI communicator to construct the interface
     *                 with. If nullptr is passed, MPI_COMM_WORLD is duplicated if
     *                 CiE was compiled with MPI support, otherwise a serial inter-
     *                 face is constructed.
     *  @details If the static MPI interface was not constructed and initialized yet,
     *           it will be after this function is called. The constructed interface
     *           is distributed if the input argument @a distributed is true, otherwise
     *           a serial dummy interface is created.
     *  @throws If the static MPI interface already exist and @a p_comm is not its communicator.
     */
    static Ref<MPI> get(Ptr<MPI_Comm> p_comm);

    /** @brief Get the static MPI interface.
     *
     *  @details If the static MPI interface was not constructed and initialized yet,
     *           it will be after this function is called. If CiE was not compiled
     *           wiht MPI support, a serial dummy MPI interface is constructed.
     *  @note If you need to control what kind of MPI interface to construct, use
     *        @ref MPISingleton::get(Ptr<MPI_Comm>).
     */
    static std::shared_ptr<MPI> share();

    /** @brief Get the static MPI interface.
     *
     *  @param p_comm: Ptr to the MPI communicator to construct the interface
     *                 with. If nullptr is passed, MPI_COMM_WORLD is duplicated if
     *                 CiE was compiled with MPI support, otherwise a serial inter-
     *                 face is constructed.
     *  @details If the static MPI interface was not constructed and initialized yet,
     *           it will be after this function is called. The constructed interface
     *           is distributed if the input argument @a distributed is true, otherwise
     *           a serial dummy interface is created.
     *  @throws If the static MPI interface already exist and @a p_comm is not its communicator.
     */
    static std::shared_ptr<MPI> share(Ptr<MPI_Comm> p_comm);

    /// @brief Kill the static MPI interface if possible.
    /// @throws If this singleton isn't the sole owner of the static MPI interface.
    static void reset();

private:
    static void ensureInitialized(bool req, Ptr<MPI_Comm> p_comm);

    static std::optional<std::shared_ptr<MPI>> _p_mpi;

    static Mutex _mutex;
}; // class MPISingleton


} // namespace cie::mpi

#endif
