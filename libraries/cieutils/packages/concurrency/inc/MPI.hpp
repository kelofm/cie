#ifndef CIE_MPI_HPP
#define CIE_MPI_HPP

// --- Utility Includes ---
#include "packages/concurrency/inc/mpi_fwd.hpp"
#include "packages/io/inc/Serializer.hpp"
#include "packages/concurrency/inc/RankID.hpp"
#include "packages/concurrency/inc/MessageTag.hpp"
#include "packages/types/inc/types.hpp"

// --- STL Includes ---
#include <memory>
#include <string>
#include <optional>


namespace cie::mpi {


class MPIImpl;


class MPI
{
private:
    struct MPIKey {};

public:
    /** @brief Construct the root MPI interface.
     *
     *  @param p_comm: Ptr to the MPI_Comm to use, or nullptr.
     *  @param key: private-access key ensuring that this constructor cannot be called by outsiders.
     *  @details If @a p_comm is not null, this constructor initializes MPI and duplicates
     *           the passed communicator.
     */
    MPI(Ptr<MPI_Comm> p_comm, MPIKey key);

    ~MPI();

    /// @name Value Return Interface
    /// @{

    template <class T>
    T receive(RankID source, MessageTag tag = 0);

    /// @}
    /// @name Ref Interface
    /// @{

    template <concepts::TriviallySerializable T>
    void send(Ref<const T> r_message,
              RankID destination,
              MessageTag tag = 0);

    template <concepts::TriviallySerializable T>
    void receive(Ref<T> r_message,
                 RankID source,
                 MessageTag tag = 0);

    template <concepts::TriviallySerializable T>
    void sendAndReceive(Ref<const T> r_send,
                        RankID sendTo,
                        Ref<T> r_receive,
                        RankID receiveFrom,
                        MessageTag tag = 0);

    template <concepts::TriviallySerializable T>
    void broadcast(Ref<T> r_message,
                   RankID source,
                   MessageTag tag = 0);

    /// @}
    /// @name Synchronization
    /// @{

    void barrier();

    /// @}
    /// @name Queries
    /// @{

    Size size() const;

    RankID getRankID() const noexcept;

    RankID getMasterRankID() const;

    /// @brief Get the static MPI interface.
    static Ref<MPI> getRoot();

    Ptr<MPI_Comm> getComm();

    bool isMaster() const;

    bool isDistributed() const noexcept;

    bool isRoot() const noexcept;

    /// @}

private:
    MPI(RightRef<std::unique_ptr<MPIImpl>> rp_impl);

    MPI(MPI&&) = default;

    MPI(const MPI&&) = delete;

    MPI& operator=(MPI&&) = delete;

    MPI& operator=(const MPI&) = delete;

private:
    friend class MPIImpl;

    friend class MPISerialImpl;

    friend class MPIDistributedImpl;

    friend class MPISingleton;

    std::unique_ptr<MPIImpl> _p_impl;
}; // class MPI


} // namespace cie::mpi

#include "packages/concurrency/impl/MPI_impl.hpp"

#endif
